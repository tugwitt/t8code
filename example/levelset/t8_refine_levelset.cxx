/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <sc_refcount.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>
#include <p4est_connectivity.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_cmesh/t8_cmesh_partition.h>
#include <t8_cmesh_readmshfile.h>
#include <t8_forest.h>
#include <t8_default_cxx.hxx>
#include <example/common/t8_example_common.h>

/* Create a cmesh from a .msh files uniform level 0
 * partitioned. */
static void
t8_refine_ls_forest_cmesh_mshfile (t8_cmesh_t cmesh, const char *vtu_prefix,
                              sc_MPI_Comm comm, int init_level, int max_level,
                              int no_vtk, double midpoint[3], double radius,
                              double bandwidth, double T,
                              double delta_t, int do_ghost, int do_balance)
{
  t8_cmesh_t          cmesh_partition;
  char                forest_vtu[BUFSIZ], cmesh_vtu[BUFSIZ];
  t8_forest_t         forest, forest_adapt, forest_partition, forest_balance;
  double              t;
  int                 partition_cmesh, r;
  const int           refine_rounds = max_level - init_level;
  int                 time_step;
  t8_levelset_sphere_data_t sphere;
  t8_example_level_set_struct_t levelset_description;

  t8_global_productionf ("Committed cmesh with"
                         " %lli global trees.\n",
                         (long long) t8_cmesh_get_num_trees (cmesh));


  /* If the input cmesh is partitioned then we use a partitioned cmesh
   * and also repartition it in each timestep (happens automatically in
   * t8_forest_commit). We have to initially start with a uniformly refined
   * cmesh in order to be able to construct the forest on it.
   * If on the other hand, the input cmesh was replicated, then we keep it
   * as replicated throughout. */
  partition_cmesh = t8_cmesh_is_partitioned (cmesh);
  if (partition_cmesh) {
    /* Set up cmesh_partition to be a repartition of cmesh. */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    /* The new cmesh is partitioned according to a uniform init_level refinement */
    t8_cmesh_set_partition_uniform (cmesh_partition, init_level);
    t8_cmesh_commit (cmesh_partition, comm);
  }
  else {
    /* Use cmesh_partition as the original replicated cmesh */
    cmesh_partition = cmesh;
  }
  /* Initialize forest and set cmesh */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh_partition, comm);
  /* Set the element scheme */
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  /* Set the initial refinement level */
  t8_forest_set_level (forest, init_level);
  /* Commit the forest */
  t8_forest_commit (forest);

  /* vtu output */
  if (!no_vtk) {
    snprintf (forest_vtu, BUFSIZ, "%s_forest_uniform", vtu_prefix);
    t8_forest_write_vtk (forest, forest_vtu);
    t8_debugf ("Wrote adapted forest\n");
  }

  /* Set the data for adapt. */
  /* Set the levelset description */
  levelset_description.L = t8_levelset_sphere; /* Refinement along a sphere */
  levelset_description.band_width = bandwidth;
  levelset_description.max_level = max_level;
  levelset_description.min_level = init_level;
  levelset_description.t = 0; /* initial time is 0 */
  levelset_description.udata = &sphere; /* The sphere around which we refine */
  /* Define the sphere */
  sphere.M[0] = midpoint[0];
  sphere.M[1] = midpoint[1];
  sphere.M[2] = midpoint[2];
  sphere.radius = radius;

  /* Start the time loop, in each time step the refinement moves
   * further through the domain */
  for (t = 0, time_step = 0; t < T; t += delta_t, time_step++) {
    /* Adapt the forest */
    for (r = 0; r < refine_rounds; r++) {
      t8_forest_init (&forest_adapt);
      t8_forest_set_adapt (forest_adapt, forest, t8_common_adapt_level_set, 0);
      /* Move the sphere with time in X-direction */
      sphere.M[0] += t;
      t8_forest_set_user_data (forest_adapt, (void *) &levelset_description);
      t8_global_productionf ("Starting levelset refinement around {%.2f,%.2f,%.2f} "
                             "with radius %.3f, band width %.3f\n", sphere.M[0],
                             sphere.M[1], sphere.M[2], sphere.radius, bandwidth);
      t8_global_productionf ("Levels: %i to %i\n", init_level, max_level);

      t8_forest_commit (forest_adapt);
      forest = forest_adapt;
      forest_adapt = NULL;
    }
      /* vtu output */
      if (!no_vtk) {
        snprintf (forest_vtu, BUFSIZ, "%s_forest_adapt_%03d", vtu_prefix,
                  time_step);
        t8_forest_write_vtk (forest, forest_vtu);
        t8_debugf ("Wrote adapted forest\n");
      }


      /* If desired, balance after last step */
      if (do_balance) {
        t8_forest_init (&forest_balance);
        t8_forest_set_balance (forest_balance, forest, 0);
        t8_forest_commit (forest_balance);
        forest = forest_balance;

        /* vtu output */
        if (!no_vtk) {
          snprintf (forest_vtu, BUFSIZ, "%s_forest_balance_%03d", vtu_prefix,
                    time_step);
          t8_forest_write_vtk (forest, forest_vtu);
          t8_debugf ("Wrote adapted forest\n");
        }
      }

      /* partition the adapted forest */
      t8_forest_init (&forest_partition);
      /* partition the adapted forest */
      t8_forest_set_partition (forest_partition, forest, 0);
      if (do_ghost) {
        t8_forest_set_ghost (forest_partition, 1, T8_GHOST_FACES);
      }
      t8_forest_commit (forest_partition);

      /* vtu output */
      if (!no_vtk) {
        snprintf (forest_vtu, BUFSIZ, "%s_forest_partition_%03d", vtu_prefix,
                  time_step);
        t8_forest_write_vtk (forest_partition, forest_vtu);
        t8_debugf ("Wrote adapted forest\n");
      }

      forest = forest_partition;

    /* Set the vtu output name */
    if (!no_vtk) {
      snprintf (forest_vtu, BUFSIZ, "%s_forest_partition_%03d", vtu_prefix,
                time_step);
      snprintf (cmesh_vtu, BUFSIZ, "%s_cmesh_partition_%03d", vtu_prefix,
                time_step);
      t8_forest_write_vtk (forest_partition, forest_vtu);
      t8_cmesh_vtk_write_file (t8_forest_get_cmesh (forest_partition),
                               cmesh_vtu, 1.0);
      t8_debugf ("Wrote partitioned forest and cmesh\n");
    }
    if (partition_cmesh) {
      /* Print runtimes and statistics of forest and cmesh partition */
      t8_cmesh_print_profile (t8_forest_get_cmesh (forest_partition));
    }
    t8_forest_print_profile (forest_partition);
    /* Set forest to the partitioned forest, so it gets adapted
     * in the next time step. */
    forest = forest_partition;
    /* TIME-LOOP ends here */
  }
  /* memory clean-up */
  t8_forest_unref (&forest_partition);
}

/* Construct a cmesh either from a .msh mesh file or from a
 * collection of cmesh files constructed with t8_cmesh_save.
 * If msh_file is NULL, the cmesh is loaded from the cmesh_file and num_files
 * must be specified. If cmesh_file is NULL, the cmesh is loaded from the .msh
 * file and mesh_dim must be specified. */
t8_cmesh_t
t8_refine_ls_forest_create_cmesh (const char *msh_file, int mesh_dim,
                             const char *cmesh_file, int num_files,
                             sc_MPI_Comm comm, int init_level, int stride)
{
  t8_cmesh_t          cmesh;
  t8_cmesh_t          cmesh_partition;
  int                 partition;

  T8_ASSERT (msh_file == NULL || cmesh_file == NULL);

  if (msh_file != NULL) {
    /* Create a cmesh from the given mesh files */
    cmesh = t8_cmesh_from_msh_file ((char *) msh_file, 1, comm, mesh_dim, 0);
    partition = 1;
  }
  else {
    T8_ASSERT (cmesh_file != NULL);
    SC_CHECK_ABORT (num_files > 0, "Must specify valid number of files.\n");
    /* Load the cmesh from the stored files and evenly distribute it
     * among all ranks */
    cmesh = t8_cmesh_load_and_distribute (cmesh_file, num_files, comm,
                                          T8_LOAD_STRIDE, stride);
    /* Partition only if more than 1 input file */
    partition = num_files > 1;
  }
  SC_CHECK_ABORT (cmesh != NULL, "Error when creating cmesh.\n");

  if (partition) {
    /* partition the cmesh uniformly */
    t8_cmesh_init (&cmesh_partition);
    t8_cmesh_set_derive (cmesh_partition, cmesh);
    t8_cmesh_set_partition_uniform (cmesh_partition, init_level);
    t8_cmesh_commit (cmesh_partition, comm);
    return cmesh_partition;
  }
  return cmesh;
}

int
main (int argc, char *argv[])
{
  int                 mpiret, mpisize;
  int                 first_argc;
  int                 level, level_diff;
  int                 help = 0, no_vtk, do_ghost, do_balance;
  int                 dim, num_files;
  int                 stride;
  double              T, delta_t, cfl;
  sc_options_t       *opt;
  t8_cmesh_t          cmesh;
  const char         *mshfileprefix, *cmeshfileprefix; /* prefix of the input files */
  const char         *vtu_prefix; /* prefix of the output files */
  double              sphere_midpoint[3]; /* Midpoint of the sphere around which we refine */
  double              sphere_radius; /* Radius of the sphere */
  double              bandwidth; /* bandwidth of the refinement region */

  /* Initialize MPI, sc, p4est and t8code */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* get mpisize */
  mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &mpisize);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_STATISTICS);

  /* Setup for command line options */
  opt = sc_options_new (argv[0]);

  sc_options_add_switch (opt, 'h', "help", &help,
                         "Display a short help message.");
  sc_options_add_switch (opt, 'o', "no-vtk", &no_vtk,
                         "Do not write vtk output.");
  sc_options_add_string (opt, 'f', "mshfile", &mshfileprefix, NULL,
                         "If specified, the cmesh is constructed from a .msh file with "
                         "the given prefix. The files must end in .msh and be "
                         "created with gmsh. If neither -f or -c are used, "
                         "a unit-square quad mesh is taken as cmesh.");
  sc_options_add_int (opt, 'd', "dim", &dim, 2,
                      "Together with -f: The dimension of the coarse mesh. 2 or 3.");
  sc_options_add_string (opt, 'c', "cmeshfile", &cmeshfileprefix, NULL,
                         "If specified, the cmesh is constructed from a collection "
                         "of cmesh files. Created with t8_cmesh_save."
                         "The number of files must then be specified with the -n "
                         "option. If neither -f or -c are used, "
                         "a unit-square quad mesh is taken as cmesh.");
  sc_options_add_int (opt, 'n', "nfiles", &num_files, -1,
                      "If the -c option is used, the number of cmesh files must "
                      "be specified as an argument here. If n=1 then the cmesh "
                      "will be replicated throughout the test.");
  sc_options_add_int (opt, 's', "stride", &stride, 16,
                      "If -c and -n are used, only every s-th MPI rank will "
                      "read a .cmesh file (file number: rank/s). Default is 16.");
  sc_options_add_int (opt, 'l', "level", &level, 0,
                      "The initial uniform "
                      "refinement level of the forest.");
  sc_options_add_int (opt, 'r', "rlevel", &level_diff, 1,
                      "The number of levels that the forest is refined "
                      "from the initial level.");
  sc_options_add_double (opt, 'x', "midpoint-x", sphere_midpoint, 0,
                         "The x coordinate of the sphere's midpoint.");
  sc_options_add_double (opt, 'y', "midpoint-y", sphere_midpoint + 1, 0,
                         "The y coordinate of the sphere's midpoint.");
  sc_options_add_double (opt, 'z', "midpoint-z", sphere_midpoint + 2, 0,
                         "The z coordinate of the sphere's midpoint.");
  sc_options_add_double (opt, 'R', "radius", &sphere_radius, 1,
                         "The radius of the sphere.");
  sc_options_add_double (opt, 'B', "bandwidth", &bandwidth, 1,
                         "The bandwidth of the refinement region.");
  sc_options_add_double (opt, 'T', "time", &T, 1,
                         "The simulated time span."
                         "We simulate the time from 0 to T");
  sc_options_add_double (opt, 'D', "delta_t", &delta_t, 0.08,
                         "The time step in each simulation step. "
                         "Deprecated, use -C instead.");
  /* CFL number. delta_t = CFL * 0.64 / 2^level */
  sc_options_add_double (opt, 'C', "cfl", &cfl, 0,
                         "The CFL number. If specified, then delta_t is set to "
                         "CFL * 0.64 / 2^level. Overwrites any other delta_t setting.");
  sc_options_add_switch (opt, 'g', "ghost", &do_ghost,
                         "Create ghost elements.");
  sc_options_add_switch (opt, 'b', "balance", &do_balance,
                         "Establish a 2:1 balance in the forest.");

  /* parse command line options */
  first_argc = sc_options_parse (t8_get_package_id (), SC_LP_DEFAULT,
                                 opt, argc, argv);
  /* check for wrong usage of arguments */
  if (first_argc < 0 || first_argc != argc || dim < 2 || dim > 3
      || stride <= 0
      || (num_files - 1) * stride >= mpisize || cfl < 0 || sphere_radius < 0
      || bandwidth < 0) {
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
    return 1;
  }
  if (help) {
    /* Display help message */
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else {
    /* Execute this part of the code if all options are correctly set */
    /* Set correct timestep, if -C was specified with a value greater 0, then
     * we overwrite the delta_t setting. We choose this to be backwards compatible to
     * call that use the delta_t option. Eventually, we will remove the delta_t option
     * completely. */
    if (cfl > 0) {
      delta_t = cfl * 0.64 / (1 << level);
    }
    t8_global_productionf ("Using delta_t = %f\n", delta_t);
    if (mshfileprefix != NULL) {
      cmesh = t8_refine_ls_forest_create_cmesh (mshfileprefix, dim, NULL, -1,
                                           sc_MPI_COMM_WORLD, level, stride);
      vtu_prefix = mshfileprefix;
    }
    else if (cmeshfileprefix != NULL) {
      cmesh = t8_refine_ls_forest_create_cmesh (NULL, -1, cmeshfileprefix,
                                           num_files, sc_MPI_COMM_WORLD,
                                           level, stride);
      vtu_prefix = cmeshfileprefix;
    }
    else {
      cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);
      vtu_prefix = "unitsquare";
    }
    t8_refine_ls_forest_cmesh_mshfile (cmesh, vtu_prefix,
                                  sc_MPI_COMM_WORLD, level,
                                  level + level_diff, no_vtk, sphere_midpoint,
                                  sphere_radius, bandwidth, T,
                                  delta_t, do_ghost, do_balance);
  }
  sc_options_destroy (opt);
  sc_finalize ();
  return 0;
}
