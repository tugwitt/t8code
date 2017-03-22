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
#include <t8_eclass.h>
#include <t8_element_cxx.hxx>
#include <t8_default_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh.h>
#include <t8_default/t8_dtri.h>
#include <p4est.h>

/* Construct a coarse mesh of two trees that are connected to each other.
 * Build a forest on top of it and perform some element_neighbor routines.
 * If hybrid is true, eclass is ignored and a hybrid quad/triangle mesh is created. */
void
t8_ghost_neighbor_test (t8_eclass_t eclass, sc_MPI_Comm comm, int hybrid)
{
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_element_t       *elem, *neigh1, *neigh2, *neigh;
  t8_scheme_cxx_t    *scheme;
  t8_eclass_scheme_c *ts, *triangle_scheme, *quad_scheme;
  int                 level = 1;
  t8_locidx_t         element_id = 1;

  if (hybrid) {
    eclass = T8_ECLASS_QUAD;
  }
  T8_ASSERT (eclass != T8_ECLASS_VERTEX);
  t8_cmesh_init (&cmesh);
  /* We construct a coarse mesh of 2 trees that are connected to each other. */
  t8_cmesh_set_tree_class (cmesh, 0, eclass);
  t8_cmesh_set_tree_class (cmesh, 1, hybrid ? T8_ECLASS_TRIANGLE : eclass);
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 1);
  t8_cmesh_commit (cmesh, comm);

  /* Construct a forest */
  t8_forest_init (&forest);
  t8_forest_set_cmesh (forest, cmesh, comm);
  t8_forest_set_level (forest, level);
  /* construct the scheme and get the eclass scheme */
  scheme = t8_scheme_new_default_cxx ();
  ts = scheme->eclass_schemes[eclass];
  t8_forest_set_scheme (forest, scheme);
  t8_forest_commit (forest);

  /* hybrid does only work with element 0 or 1 */
  T8_ASSERT (!hybrid
             || (element_id == 0 || element_id == 1 || element_id == 3)
             && level == 1);
  /* Get the data of element number element_id */
  elem = t8_forest_get_element (forest, element_id);
  T8_ASSERT (elem != NULL);
  /* allocate memory for the face neighbor */
  if (hybrid) {
    triangle_scheme = scheme->eclass_schemes[T8_ECLASS_TRIANGLE];
    ts = quad_scheme = scheme->eclass_schemes[T8_ECLASS_QUAD];
    triangle_scheme->t8_element_new (1, &neigh2);
  }
  ts->t8_element_new (1, &neigh1);
  {
    /* Iterate over all faces and create the face neighbor */
    int                 i, ret;
    int                 anchor_node[3];

    t8_debugf ("root len = %i\n", ts->t8_element_root_len (elem));
    for (i = 0; i < t8_eclass_num_faces[eclass]; i++) {
      neigh = neigh1;
      if (hybrid && i == 1 && element_id == 1) {
        /* face neighbor is in triangle tree */
        ts = triangle_scheme;
        neigh = neigh2;
      }
      else if (hybrid) {
        ts = quad_scheme;
      }
      ret = t8_forest_element_face_neighbor (forest, 0, elem, neigh, i);
      if (ret != -1) {
        ts->t8_element_anchor (neigh, anchor_node);
        t8_debugf
          ("neighbor of %i across face %i (in tree %i): (%i,%i,%i,%i)\n",
           element_id, i,
           ret, ts->t8_element_level (neigh), anchor_node[0],
           anchor_node[1], t8_eclass_to_dimension[eclass] > 2 ? anchor_node[2]
           : -1);
      }
    }
  }
  if (hybrid) {
    triangle_scheme->t8_element_destroy (1, &neigh2);
    ts = quad_scheme;
  }
  ts->t8_element_destroy (1, &neigh1);
  t8_forest_unref (&forest);
}

int
main (int argc, char **argv)
{
  int                 mpiret;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  t8_global_productionf ("Testing neighbors for triangle\n");
  t8_ghost_neighbor_test (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD, 0);
  t8_global_productionf ("Testing neighbors for quad\n");
  t8_ghost_neighbor_test (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0);
  t8_global_productionf ("Testing neighbors for tet\n");
  t8_ghost_neighbor_test (T8_ECLASS_TET, sc_MPI_COMM_WORLD, 0);
  t8_global_productionf ("Testing neighbors for hex\n");
  t8_ghost_neighbor_test (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0);
  t8_global_productionf ("Testing neighbors for hybrid_mesh\n");
  t8_ghost_neighbor_test (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 1);

  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}