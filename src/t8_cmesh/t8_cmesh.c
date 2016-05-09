/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

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

#include <t8_cmesh.h>
#include <t8_cmesh_vtk.h>
#include <t8_refcount.h>
#include <t8_shmem.h>
#ifdef T8_WITH_METIS
#include <metis.h>
#endif
#include "t8_cmesh_trees.h"

/** \file t8_cmesh.c
 *
 * TODO: document this file
 */

int
t8_cmesh_is_initialized (t8_cmesh_t cmesh)
{
  if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) &&
        !cmesh->committed)) {
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* TODO: check conditions that must always hold after init and before commit */
  if (0) {
    return 0;
  }
#endif

  return 1;
}

int
t8_cmesh_is_committed (t8_cmesh_t cmesh)
{
  if (!(cmesh != NULL && t8_refcount_is_active (&cmesh->rc) &&
        cmesh->committed)) {
    return 0;
  }

#ifdef T8_ENABLE_DEBUG
  /* TODO: check more conditions that must always hold after commit */
  if ((!t8_cmesh_trees_is_face_consistend (cmesh, cmesh->trees)) || 0) {
    return 0;
  }
#endif

  return 1;
}

/* *INDENT-OFF* */
static t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_locidx_t tree_id);
/* *INDENT-ON* */

#if 0
/* Compute a hash value for a ghost tree. */
/* deprecated */
static unsigned
t8_cmesh_ghost_hash_fn (const void *ghost, const void *data)
{
  t8_cmesh_t          cmesh;
  t8_cghost_t         G;

  T8_ASSERT (data != NULL);
  cmesh = (t8_cmesh_t) data;
  T8_ASSERT (cmesh->num_ghosts > 0);
  T8_ASSERT (cmesh->set_partition);
  T8_ASSERT (cmesh->num_local_trees > 0);

  G = (t8_cghost_t) ghost;
  /* TODO: is this a reasonable hash value? */
  return G->treeid % cmesh->num_ghosts;
}

static int
t8_cmesh_ghost_equal_fn (const void *ghost1, const void *ghost2,
                         const void *data)
{
  t8_cghost_t         G1, G2;

  G1 = (t8_cghost_t) ghost1;
  G2 = (t8_cghost_t) ghost2;

  return G1->treeid == G2->treeid;
}
#endif

/* Check whether a given communicator assigns the same rank and mpisize
 * as stored in a given cmesh. */
int
t8_cmesh_comm_is_valid (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 mpiret, mpisize, mpirank;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  if (mpisize != cmesh->mpisize || mpirank != cmesh->mpirank) {
    return 0;
  }
  return 1;
}

void
t8_cmesh_init (t8_cmesh_t * pcmesh)
{
  t8_cmesh_t          cmesh;
  T8_ASSERT (pcmesh != NULL);

  cmesh = *pcmesh = T8_ALLOC_ZERO (t8_cmesh_struct_t, 1);
  t8_refcount_init (&cmesh->rc);

  /* sensible (hard error) defaults */
  cmesh->set_refine_level = 0;  /*< sensible default TODO document */
  cmesh->set_partition_level = -1;
  cmesh->dimension = -1;        /*< ok; force user to select dimension */
  cmesh->mpirank = -1;
  cmesh->mpisize = -1;
  cmesh->face_knowledge = 3;    /*< sensible default TODO document */
  t8_stash_init (&cmesh->stash);

  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
}

#if 0
/* This function is not part of the interface. The number of trees is always clear
 * from the number of calls to t8_cmesh_set_tree_class.
 * It is set in t8_cmesh_commit */
/* TODO: rename num_trees to global_num_trees or num_gtrees etc.
 *       to always distinguish between local and global.
 *       Do this everywhere in the code.
 */
static void
t8_cmesh_set_num_trees (t8_cmesh_t cmesh, t8_gloidx_t num_trees)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  /* If the cmesh is entered as a partitioned cmesh,
   * this function sets the local number of trees;
   * (TODO ^^^ would require locidx -- better provide two arguments)
   * the global number then must have been set in cmesh_set_partition.
   * Otherwise the global number of trees is set here.
   * TODO: make this function behave consistently independent on prior
   *       calls to set_partition.
   *       We want the user to be free in the sequence of calls
   *       as much as possible.
   */
  if (cmesh->set_partition) {
    /* num_trees == 0 is allowed */
    T8_ASSERT (cmesh->num_trees > 0);
    T8_ASSERT (cmesh->num_local_trees == 0);
    cmesh->num_local_trees = num_trees;
  }
  else {
    /* num_trees == 0 is allowed */
    T8_ASSERT (cmesh->num_trees >= 0);
    cmesh->num_trees = cmesh->num_local_trees = num_trees;
  }
  /* As soon as we know the number of trees, we allocate
   * the ctree array.
   * TODO?
   */
}
#endif

void
t8_cmesh_set_derive (t8_cmesh_t cmesh, t8_cmesh_t set_from)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (set_from == NULL || t8_cmesh_is_committed (set_from));

  /* TODO: If a previuously set cmesh is overwritten than the user is
   * responsible for unreffing it. Otherwise we have to give an mpi communicator
   * as parameter to this function. */
  cmesh->set_from = set_from;
  if (set_from != NULL) {
    t8_cmesh_ref (set_from);
  }
}

#if 0
/* TODO: deprecated, remove */
void
t8_cmesh_set_partition (t8_cmesh_t cmesh, int set_partition,
                        int set_face_knowledge,
                        t8_gloidx_t first_local_tree,
                        t8_gloidx_t last_local_tree,
                        t8_gloidx_t * tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (0 <= set_face_knowledge && set_face_knowledge <= 3);
  /* TODO: allow -1 for set_face_knowledge to keep it unchanged?
   *      update: unchanged from what? face_knowledge is only important for the
   * information on the stash. When the cmesh is derived there is no
   * stash. A committed cmesh has always face_knowledge 3. */

  /* TODO: Careful with tese assumptions; allow the user maximum flexibility */
#if 0
  T8_ASSERT (cmesh->num_trees == 0);
  T8_ASSERT (cmesh->num_local_trees == 0);
  T8_ASSERT (cmesh->first_tree == 0);
#endif

  /* set cmesh->set_partition to 0 or 1 (no; we always treat nonzero as true) */
  cmesh->set_partition = set_partition;
  /* TODO: this is how to query boolean variables */
  if (set_partition) {
    cmesh->first_tree = first_local_tree;
    cmesh->num_local_trees = last_local_tree - first_local_tree + 1;
    /* Since num_local_trees is a locidx we have to check whether we did create an
     * overflow in the previous computation */
    T8_ASSERT (cmesh->num_local_trees ==
               last_local_tree - first_local_tree + 1);
    cmesh->face_knowledge = set_face_knowledge;
    /* Right now no other face_knowledge is supported */
    SC_CHECK_ABORTF (set_face_knowledge == 3, "Level %i of face knowledge"
                     "is not supported.\n", set_face_knowledge);
    cmesh->tree_offsets = tree_offsets;
  }
}
#endif

void
t8_cmesh_set_partition_range (t8_cmesh_t cmesh, int set_face_knowledge,
                              t8_gloidx_t first_local_tree,
                              t8_gloidx_t last_local_tree)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  SC_CHECK_ABORT (set_face_knowledge == 1 || set_face_knowledge == 3,
                  "Face knowledge other than three is not implemented yet.");
  cmesh->face_knowledge = set_face_knowledge;
  cmesh->first_tree = first_local_tree;
  cmesh->num_local_trees = last_local_tree - first_local_tree + 1;
  cmesh->set_partition = 1;
  /* Overwrite previous partition settings */
  cmesh->tree_offsets = NULL;   /* TODO: Do we need to free memory? */
  cmesh->set_partition_level = -1;
}

void
t8_cmesh_set_partition_offsets (t8_cmesh_t cmesh, t8_gloidx_t * tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));

  cmesh->tree_offsets = tree_offsets;
  cmesh->set_partition = 1;
  if (tree_offsets != NULL) {
    /* We overwrite any previously partition settings */
    cmesh->first_tree = -1;
    cmesh->num_local_trees = -1;
    cmesh->set_partition_level = -1;
  }
}

void
t8_cmesh_set_partition_uniform (t8_cmesh_t cmesh, int element_level)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (element_level >= -1);

  cmesh->set_partition = 1;
  cmesh->set_partition_level = element_level;
  if (element_level >= 0) {
    /* We overwrite any previous partition settings */
    cmesh->first_tree = -1;
    cmesh->num_local_trees = -1;
    cmesh->tree_offsets = NULL; /* TODO: Do we need to free memory? */
  }
}

#if 0
/* No longer needed */
void
t8_cmesh_set_partition_from (t8_cmesh_t cmesh, const t8_cmesh_t cmesh_from,
                             int level, t8_gloidx_t * tree_offsets)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (t8_cmesh_is_committed (cmesh_from));
  T8_ASSERT (cmesh_from->set_partition);

  cmesh->set_from = cmesh_from;
  cmesh->set_partition = 1;
  cmesh->face_knowledge = cmesh_from->face_knowledge;
  if (level >= 0) {
    cmesh->set_partition_level = level;
  }
  else {
    cmesh->tree_offsets = tree_offsets;
  }
  cmesh->from_method |= T8_CMESH_PARTITION;
}
#endif

void
t8_cmesh_set_refine (t8_cmesh_t cmesh, int level)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (level >= 0);

  cmesh->set_refine_level = level;
}

t8_gloidx_t
t8_cmesh_get_first_treeid (t8_cmesh_t cmesh)
{
  return cmesh->first_tree;
}

/* Return a pointer to the ctree of a given global tree_id. */
/* TODO: should get a gloidx?
 *       place after commit */
static              t8_ctree_t
t8_cmesh_get_tree (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (0 <= tree_id && tree_id < cmesh->num_local_trees);
  T8_ASSERT (cmesh->committed);

  return t8_cmesh_trees_get_tree (cmesh->trees, tree_id);
}

/* Returns the first local tree.
 * Returns NULL if there are no local trees. */
/* TODO: hide */
t8_ctree_t
t8_cmesh_get_first_tree (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_local_trees > 0 ? t8_cmesh_get_tree (cmesh, 0) : NULL;
}

/* returns the next local tree in the cmesh (by treeid)
 * after a given tree.
 * The given tree must be a valid and owned tree.
 * If the given tree is the last local tree, NULL is returned */
/* TODO: hide */
t8_ctree_t
t8_cmesh_get_next_tree (t8_cmesh_t cmesh, t8_ctree_t tree)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (tree != NULL);
  T8_ASSERT (0 <= tree->treeid && tree->treeid < cmesh->num_local_trees);
  T8_ASSERT (cmesh->committed);
  return tree->treeid <
    cmesh->num_local_trees -
    1 ? t8_cmesh_get_tree (cmesh, tree->treeid + 1) : NULL;
}

void
t8_cmesh_set_attribute (t8_cmesh_t cmesh, t8_gloidx_t tree_id, int package_id,
                        int key, void *data, size_t data_size,
                        int data_persists)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, tree_id, package_id, key, data_size,
                          data, data_persists);
}

void               *
t8_cmesh_get_attribute (t8_cmesh_t cmesh, int package_id, int key,
                        t8_locidx_t tree_id)
{
  T8_ASSERT (cmesh->committed);
  return t8_cmesh_trees_get_attribute (cmesh->trees, tree_id, package_id,
                                       key);
}

#if 0
/* Check whether a given tree_id belongs to a tree in the cmesh.
 * If partitioned only local trees are allowed.
 */
static int
t8_cmesh_tree_id_is_owned (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  T8_ASSERT (cmesh->committed);
  if (cmesh->set_partition) {
    return cmesh->first_tree <= tree_id
      && tree_id < cmesh->first_tree + cmesh->num_local_trees;
  }
  else {
    return 0 <= tree_id && tree_id < cmesh->num_trees;
  }
}

#endif

#if 0
/* Given a tree_id return the index of the specified tree in
 * cmesh's tree array
 */
static              t8_topidx_t
t8_cmesh_tree_index (t8_cmesh_t cmesh, t8_topidx_t tree_id)
{
  return cmesh->set_partition ? tree_id - cmesh->first_tree : tree_id;
}
#endif

void
t8_cmesh_set_tree_class (t8_cmesh_t cmesh, t8_gloidx_t gtree_id,
                         t8_eclass_t tree_class)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (gtree_id >= 0);

  /* If we insert the first tree, set the dimension of the cmesh
   * to this tree's dimension. Otherwise check whether the dimension
   * of the tree to be inserted equals the dimension of the cmesh. */
  if (cmesh->dimension == -1) {
    cmesh->dimension = t8_eclass_to_dimension[tree_class];
  }
  else {
    /* TODO: This makes it illegal to set a tree to i.e. quad and change it
     *       to hex later. Even if we replace all trees with another dimension.
     *       We could move this check to commit. */
    T8_ASSERT (t8_eclass_to_dimension[tree_class] == cmesh->dimension);
  }

  t8_stash_add_class (cmesh->stash, gtree_id, tree_class);
#ifdef T8_ENABLE_DEBUG
  cmesh->inserted_trees++;
#endif
}

void
t8_cmesh_set_tree_vertices (t8_cmesh_t cmesh, t8_topidx_t tree_id,
                            int package_id, int key,
                            double *vertices, t8_topidx_t num_vertices)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (vertices != NULL);
  T8_ASSERT (!cmesh->committed);

  t8_stash_add_attribute (cmesh->stash, tree_id, package_id, key,
                          3 * num_vertices * sizeof (double),
                          (void *) vertices, 1);
}

void
t8_cmesh_set_join (t8_cmesh_t cmesh, t8_gloidx_t gtree1, t8_gloidx_t gtree2,
                   int face1, int face2, int orientation)
{
  T8_ASSERT (0 <= orientation);

  t8_stash_add_facejoin (cmesh->stash, gtree1, gtree2, face1, face2,
                         orientation);
}

/* returns true if cmesh_a equals cmesh_b */
int
t8_cmesh_is_equal (t8_cmesh_t cmesh_a, t8_cmesh_t cmesh_b)
/* TODO: rewrite */
{
  int                 is_equal;
  T8_ASSERT (cmesh_a != NULL && cmesh_b != NULL);

  if (cmesh_a == cmesh_b) {
    return 1;
  }
  /* check entries that are numbers */
  is_equal = cmesh_a->committed != cmesh_b->committed || cmesh_a->dimension !=
    cmesh_b->dimension ||
    cmesh_a->set_partition != cmesh_b->set_partition ||
    cmesh_a->mpirank != cmesh_b->mpirank ||
    cmesh_a->mpisize != cmesh_b->mpisize ||
    cmesh_a->num_trees != cmesh_b->num_trees ||
    cmesh_a->num_local_trees != cmesh_b->num_local_trees ||
    cmesh_a->num_ghosts != cmesh_b->num_ghosts ||
    cmesh_a->first_tree != cmesh_b->first_tree;
#ifdef T8_ENABLE_DEBUG
  is_equal = is_equal || cmesh_a->inserted_trees != cmesh_b->inserted_trees ||
    cmesh_a->inserted_ghosts != cmesh_b->inserted_ghosts;
#endif
  if (is_equal != 0) {
    return 0;
  }
  /* check arrays */
  is_equal = memcmp (cmesh_a->num_trees_per_eclass,
                     cmesh_b->num_trees_per_eclass,
                     T8_ECLASS_COUNT * sizeof (t8_topidx_t));

  /* check tree_offsets */
  if (cmesh_a->tree_offsets != NULL) {
    if (cmesh_b->tree_offsets == NULL) {
      return 0;
    }
    else {
      is_equal = is_equal || memcmp (cmesh_a->tree_offsets,
                                     cmesh_b->tree_offsets,
                                     (cmesh_a->mpisize + 1)
                                     * sizeof (t8_gloidx_t));
    }
  }
  if (is_equal != 0) {
    return 0;
  }
  /* check trees */
  if (cmesh_a->committed &&
      !t8_cmesh_trees_is_equal (cmesh_a, cmesh_a->trees, cmesh_b->trees)) {
    /* if we have committed check tree arrays */
    return 0;
  }
  else {
    if (!cmesh_a->committed &&
        !t8_stash_is_equal (cmesh_a->stash, cmesh_b->stash)) {
      /* if we have not committed check stash arrays */
      return 0;
    }
  }
  return 1;
}

#if 0
/* broadcast the tree attributes of a cmesh on root to all processors */
/* TODO: can we optimize it by just sending the memory of the mempools? */
static void
t8_cmesh_bcast_attributes (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int                 mpirank, mpisize, mpiret;
  int                 has_attr;
  t8_ctree_t          tree;

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  for (tree = t8_cmesh_get_first_tree (cmesh_in); tree != NULL;
       tree = t8_cmesh_get_next_tree (cmesh_in, tree)) {
    if (mpirank == root && tree->attribute != NULL) {
      has_attr = 1;
    }
    else {
      has_attr = 0;
    }
    mpiret = sc_MPI_Bcast (&has_attr, 1, sc_MPI_INT, root, comm);
    SC_CHECK_MPI (mpiret);
    if (has_attr) {
      if (mpirank != root) {
        tree->attribute =
          sc_mempool_alloc (cmesh_in->tree_attributes_mem[tree->eclass]);
      }
      mpiret = sc_MPI_Bcast (tree->attribute,
                             t8_cmesh_get_attribute_size (cmesh_in,
                                                          tree->eclass),
                             sc_MPI_BYTE, root, comm);
      SC_CHECK_MPI (mpiret);
    }
  }
}
#endif

t8_cmesh_t
t8_cmesh_bcast (t8_cmesh_t cmesh_in, int root, sc_MPI_Comm comm)
{
  int                 mpirank, mpisize, mpiret;
  int                 iclass;

  struct
  {
    int                 dimension;
    t8_topidx_t         num_trees;
    t8_topidx_t         num_trees_per_eclass[T8_ECLASS_COUNT];
    size_t              stash_elem_counts[3];
#ifdef T8_ENABLE_DEBUG
    t8_topidx_t         inserted_trees;
    sc_MPI_Comm         comm;
#endif
  } dimensions;

  /* TODO: BUG: running with two processes and a cmesh of one T8_ECLASS_LINE,
   *       the on both processes the face_neigbors and vertices arrays of
   *       the single tree point to the same physical memory.
   *       (face_neighbors on both processes are equal and vertices on both
   *        processes are equal)
   */
  /* TODO: Send the tree's vertices */

  /* TODO: rewrite */

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);

  T8_ASSERT (0 <= root && root < mpisize);
  T8_ASSERT (mpirank == root || cmesh_in == NULL);
  T8_ASSERT (mpirank != root || cmesh_in != NULL);
  T8_ASSERT (mpirank != root || cmesh_in->set_partition == 0);
  /* The cmesh on the calling process must not be owned by something
   * else. */
  /* TODO: would it be useful to allow bcast even if the cmesh is referenced?
   * But then the bcasted version on other procs would have a different refcount
   * than the cmesh on the root */
  T8_ASSERT (mpirank != root || cmesh_in->rc.refcount == 1);

  /* At first we broadcast all meta information. */
  if (mpirank == root) {
    dimensions.dimension = cmesh_in->dimension;
    dimensions.num_trees = cmesh_in->num_trees;
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      dimensions.num_trees_per_eclass[iclass] =
        cmesh_in->num_trees_per_eclass[iclass];
    }
    dimensions.stash_elem_counts[0] = cmesh_in->stash->attributes.elem_count;
    dimensions.stash_elem_counts[1] = cmesh_in->stash->classes.elem_count;
    dimensions.stash_elem_counts[2] = cmesh_in->stash->joinfaces.elem_count;
#ifdef T8_ENABLE_DEBUG
    dimensions.inserted_trees = cmesh_in->inserted_trees;
#endif
  }
  /* TODO: we could optimize this by using IBcast */
  mpiret = sc_MPI_Bcast (&dimensions, sizeof (dimensions), sc_MPI_BYTE, root,
                         comm);
  SC_CHECK_MPI (mpiret);

  /* If not root store information in new cmesh and allocate memory for arrays. */
  if (mpirank != root) {
    t8_cmesh_init (&cmesh_in);
    cmesh_in->dimension = dimensions.dimension;
    cmesh_in->num_trees = dimensions.num_trees;
    for (iclass = 0; iclass < T8_ECLASS_COUNT; iclass++) {
      cmesh_in->num_trees_per_eclass[iclass] =
        dimensions.num_trees_per_eclass[iclass];
    }
#ifdef T8_ENABLE_DEBUG
    cmesh_in->inserted_trees = dimensions.inserted_trees;
    T8_ASSERT (dimensions.comm == comm);
#endif
  }
  /* broadcast all the stashed information about trees/neighbors/attributes */
  t8_stash_bcast (cmesh_in->stash, root, comm, dimensions.stash_elem_counts);
  return cmesh_in;
}

#ifdef T8_WITH_METIS
void
t8_cmesh_reorder (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  int                 mpisize, mpiret;
  idx_t               idx_mpisize;
  idx_t               ncon = 1, elemens;
  idx_t               volume, *partition, ipart, newpart;
  int                 num_faces, iface, count_face;
  idx_t              *xadj, *adjncy;
  int                 success;
  t8_locidx_t        *new_number, itree, *tree_per_part_off, *tree_per_part;
  t8_locidx_t        *face_neighbor;
  t8_locidx_t         neigh_id;
  t8_ctree_t          tree;

  /* cmesh must be commited and not partitioned */
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (!cmesh->set_partition);

  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  idx_mpisize = mpisize;
  SC_CHECK_MPI (mpiret);

  elemens = cmesh->num_trees;
  T8_ASSERT ((t8_locidx_t) elemens == cmesh->num_trees);

  /* Count the number of tree-to-tree connections via a face */
  num_faces = 0;
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor,
                                        NULL);
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      if (face_neighbor[iface] >= 0)
        num_faces++;
    }
  }

  /* xadj and adjncy store the face-connections in a CSR format
   * xadj[treeid] = offset of the tree in adjncy
   * adjncy[xadj[treeid]]...adjncy[xadj[treeid]-1] are the trees with which
   * the tree has a face connection */
  xadj = T8_ALLOC_ZERO (idx_t, elemens + 1);
  adjncy = T8_ALLOC (idx_t, num_faces);

  /* fill xadj and adjncy arrays */
  for (itree = 0, count_face = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_get_tree (cmesh, itree);
    xadj[itree + 1] = xadj[itree];
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      if (face_neighbor[iface] >= 0) {
        adjncy[count_face++] = face_neighbor[iface];
        xadj[itree + 1]++;
      }
    }
  }

  /* partutuib stores the new partitino number for each element */
  partition = T8_ALLOC (idx_t, elemens);
  /* partition the elements in mpisize many partitions */
  success =
    METIS_PartGraphRecursive (&elemens, &ncon, xadj, adjncy, NULL, NULL, NULL,
                              &idx_mpisize, NULL, NULL, NULL, &volume,
                              partition);
  T8_ASSERT (success == METIS_OK);
  /* memory to store the new treeid of a tree */
  new_number = T8_ALLOC (t8_locidx_t, cmesh->num_trees);
  /* Store the number of trees pinter partition */
  tree_per_part = T8_ALLOC_ZERO (t8_locidx_t, mpisize);
  /* Store the treeid offset of each partition. */
  tree_per_part_off = T8_ALLOC_ZERO (t8_locidx_t, mpisize + 1);
  tree_per_part_off[0] = 0;
  /* compute tree_per_part and prepare tree_per_part_off */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree_per_part[partition[itree]]++;
    tree_per_part_off[partition[itree] + 1]++;
  }
  /* compute tree_per_part_off */
  for (ipart = 1; ipart <= mpisize; ipart++) {
    tree_per_part_off[ipart] += tree_per_part_off[ipart - 1];
  }
  /* Compute for each tree its new treeid */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    newpart = partition[itree];
    T8_ASSERT (tree_per_part[newpart] > 0);
    new_number[itree] =
      tree_per_part_off[newpart + 1] - tree_per_part[newpart];
    tree_per_part[newpart]--;
  }
  /* Set for each tree its new treeid and the new ids of its neighbors */
  for (itree = 0; itree < cmesh->num_trees; itree++) {
    tree = t8_cmesh_trees_get_tree_ext (cmesh->trees, itree, &face_neighbor,
                                        NULL);
    tree->treeid = new_number[itree];
    for (iface = 0; iface < t8_eclass_num_faces[tree->eclass]; iface++) {
      neigh_id = face_neighbor[iface];
      if (neigh_id >= 0) {
        face_neighbor[iface] = new_number[neigh_id];
      }
    }
  }
  T8_FREE (partition);
  T8_FREE (xadj);
  T8_FREE (adjncy);
  T8_FREE (new_number);
  T8_FREE (tree_per_part);
  T8_FREE (tree_per_part_off);
}
#endif

static void
t8_cmesh_free_treecount (t8_cmesh_t cmesh, sc_MPI_Comm comm)
{
  SC_SHMEM_FREE (cmesh->tree_offsets, comm);
}

t8_gloidx_t
t8_cmesh_get_num_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_trees;
}

t8_locidx_t
t8_cmesh_get_num_local_trees (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  return cmesh->num_local_trees;
}

t8_eclass_t
t8_cmesh_get_tree_class (t8_cmesh_t cmesh, t8_locidx_t tree_id)
{
  t8_ctree_t          tree;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  tree = t8_cmesh_get_tree (cmesh, tree_id);
  return tree->eclass;
}

t8_eclass_t
t8_cmesh_get_ghost_class (t8_cmesh_t cmesh, t8_locidx_t ghost_id)
{
  t8_cghost_t         ghost;

  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (0 <= ghost_id && ghost_id < cmesh->num_ghosts);

  ghost = t8_cmesh_trees_get_ghost (cmesh->trees, ghost_id);
  return ghost->eclass;
}

t8_gloidx_t
t8_cmesh_get_global_id (t8_cmesh_t cmesh, t8_locidx_t local_id)
{
  T8_ASSERT (0 <= local_id && local_id <
             cmesh->num_ghosts + cmesh->num_local_trees);
  if (local_id < cmesh->num_local_trees) {
    return local_id + cmesh->first_tree;
  }
  else {
    return t8_cmesh_trees_get_ghost (cmesh->trees,
                                     local_id -
                                     cmesh->num_local_trees)->treeid;
  }
}

void
t8_cmesh_uniform_bounds (t8_cmesh_t cmesh, int level,
                         t8_gloidx_t * first_local_tree,
                         t8_gloidx_t * child_in_tree_begin,
                         t8_gloidx_t * last_local_tree,
                         t8_gloidx_t * child_in_tree_end,
                         int8_t * first_tree_shared)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);
  T8_ASSERT (level >= 0);

  *first_local_tree = 0;
  if (child_in_tree_begin != NULL) {
    *child_in_tree_begin = 0;
  }
  *last_local_tree = 0;
  if (child_in_tree_end != NULL) {
    *child_in_tree_end = 0;
  }

  if (cmesh->num_trees_per_eclass[T8_ECLASS_PYRAMID] == 0) {
    t8_gloidx_t         global_num_children;
    t8_gloidx_t         first_global_child;
    t8_gloidx_t         last_global_child;
    t8_gloidx_t         children_per_tree;
    t8_gloidx_t         prev_last_tree = -1;
    const t8_gloidx_t   one = 1;

    children_per_tree = one << cmesh->dimension * level;
    global_num_children = cmesh->num_trees * children_per_tree;

    if (cmesh->mpirank == 0) {
      first_global_child = 0;
      if (child_in_tree_begin != NULL) {
        *child_in_tree_begin = 0;
      }
    }
    else {
      /* The first global children of processor p
       * with P total processor is (the biggest int smaller than)
       * (total_num_children * p) / P
       * We cast to long double and double first to prevent integer overflow.
       */
      first_global_child =
        ((long double) global_num_children *
         cmesh->mpirank) / (double) cmesh->mpisize;
    }
    if (cmesh->mpirank != cmesh->mpisize - 1) {
      last_global_child =
        ((long double) global_num_children *
         (cmesh->mpirank + 1)) / (double) cmesh->mpisize;
    }
    else {
      last_global_child = global_num_children;
    }
    T8_ASSERT (0 <= first_global_child
               && first_global_child <= global_num_children);
    T8_ASSERT (0 <= last_global_child
               && last_global_child <= global_num_children);
    *first_local_tree = first_global_child / children_per_tree;
    if (child_in_tree_begin != NULL) {
      *child_in_tree_begin =
        first_global_child - *first_local_tree * children_per_tree;
    }
    /* TODO: Just fixed this line from last_global_child -1 / cpt
     *       Why did we not notice this error before?
     *       Changed it back*/
    *last_local_tree = (last_global_child - 1) / children_per_tree;
    if (first_tree_shared != NULL) {
      prev_last_tree = (first_global_child - 1) / children_per_tree;
      T8_ASSERT (cmesh->mpirank > 0 || prev_last_tree <= 0);
      if (cmesh->mpirank > 0 && prev_last_tree == *first_local_tree &&
          first_global_child != last_global_child && last_global_child >= 0) {
        /* We exclude empty partitions here, by def their first_tree_shared flag is zero */
        /* We also exclude that the previous partition was empty at the beginning of the
         * partitions array */
        /* TODO: If empty partitions in the middle can occur then we have to think this over */
        *first_tree_shared = 1;
      }
      else {
        *first_tree_shared = 0;
      }
    }
    if (child_in_tree_end != NULL) {
      if (*last_local_tree > 0) {
        *child_in_tree_end =
          last_global_child - *last_local_tree * children_per_tree;
      }
      else {
        *child_in_tree_end = last_global_child;
      }
    }
  }
  else {
    SC_ABORT ("Partition does not support pyramidal elements yet.");
  }
}

static void
t8_cmesh_reset (t8_cmesh_t * pcmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->rc.refcount == 0);

  /* free tree_offset */
  if (cmesh->tree_offsets != NULL) {
    T8_ASSERT (t8_cmesh_comm_is_valid (cmesh, comm));
    t8_cmesh_free_treecount (cmesh, comm);
  }
  /*TODO: write this */
  if (!cmesh->committed) {
    t8_stash_destroy (&cmesh->stash);
  }
  else {
    if (cmesh->trees != NULL) {
      t8_cmesh_trees_destroy (&cmesh->trees);
    }
  }
  if (cmesh->set_from != NULL) {
    /* We have taken ownership of cmesh_from */
    t8_cmesh_unref (&cmesh->set_from, comm);
  }

  T8_FREE (cmesh);

  *pcmesh = NULL;
}

void
t8_cmesh_ref (t8_cmesh_t cmesh)
{
  T8_ASSERT (cmesh != NULL);
  t8_refcount_ref (&cmesh->rc);
}

void
t8_cmesh_unref (t8_cmesh_t * pcmesh, sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  T8_ASSERT (pcmesh != NULL);
  cmesh = *pcmesh;
  T8_ASSERT (cmesh != NULL);

  if (t8_refcount_unref (&cmesh->rc)) {
    t8_cmesh_reset (pcmesh, comm);
  }
}

void
t8_cmesh_destroy (t8_cmesh_t * pcmesh, sc_MPI_Comm comm)
{
  T8_ASSERT (pcmesh != NULL && *pcmesh != NULL &&
             t8_refcount_is_last (&(*pcmesh)->rc));

  t8_cmesh_unref (pcmesh, comm);
  T8_ASSERT (*pcmesh == NULL);
}

/* TODO: In p4est a tree edge is joined with itself to denote a domain boundary.
 *       Will we do it the same in t8code? This is not yet decided, however the
 *       function below stores these neighbourhood information in the cmesh. */
/* TODO: Eventually we may directly partition the mesh here */
static              t8_cmesh_t
t8_cmesh_new_from_p4est_ext (void *conn, int dim, sc_MPI_Comm comm,
                             int do_dup, int set_partition)
{
#define _T8_CMESH_P48_CONN(_ENTRY) \
  (dim == 2 ? ((p4est_connectivity_t *) conn)->_ENTRY \
            : ((p8est_connectivity_t *) conn)->_ENTRY)
  t8_cmesh_t          cmesh;
  t8_topidx_t         itree;
  p4est_topidx_t      treevertex;
  double              vertices[24];     /* Only 4 * 3 = 12 used in 2d */
  int                 num_tvertices;
  int                 num_faces;
  int                 ivertex, iface;
  int8_t              ttf;
  p4est_topidx_t      ttt;

  T8_ASSERT (dim == 2 || dim == 3);
  T8_ASSERT (dim == 3 ||
             p4est_connectivity_is_valid ((p4est_connectivity_t *) (conn)));
  T8_ASSERT (dim == 2 ||
             p8est_connectivity_is_valid ((p8est_connectivity_t *) (conn)));
  num_tvertices = 1 << dim;     /*vertices per tree. 4 if dim = 2 and 8 if dim = 3. */
  num_faces = dim == 2 ? 4 : 6;
  /* basic setup */
  t8_cmesh_init (&cmesh);
  /* Add each tree to cmesh and get vertex information for each tree */
  for (itree = 0; itree < _T8_CMESH_P48_CONN (num_trees); itree++) {    /* loop over each tree */
    t8_cmesh_set_tree_class (cmesh, itree,
                             dim == 2 ? T8_ECLASS_QUAD : T8_ECLASS_HEX);
    for (ivertex = 0; ivertex < num_tvertices; ivertex++) {     /* loop over each tree corner */
      treevertex =
        _T8_CMESH_P48_CONN (tree_to_vertex[num_tvertices * itree + ivertex]);
      vertices[3 * ivertex] = _T8_CMESH_P48_CONN (vertices[3 * treevertex]);
      vertices[3 * ivertex + 1] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 1]);
      vertices[3 * ivertex + 2] =
        _T8_CMESH_P48_CONN (vertices[3 * treevertex + 2]);
    }
    t8_cmesh_set_tree_vertices (cmesh, itree, t8_get_package_id (), 0,
                                vertices, num_tvertices);
  }
  /* get face neighbor information from conn and join faces in cmesh */
  for (itree = 0; itree < _T8_CMESH_P48_CONN (num_trees); itree++) {    /* loop over each tree */
    for (iface = 0; iface < num_faces; iface++) {       /* loop over each face */
      ttf = _T8_CMESH_P48_CONN (tree_to_face[num_faces * itree + iface]);
      ttt = _T8_CMESH_P48_CONN (tree_to_tree[num_faces * itree + iface]);
      /* insert the face only if we did not insert it before */
      if (itree < ttt || (itree == ttt && iface < ttf % num_faces)) {
        t8_cmesh_set_join (cmesh, itree, ttt, iface, ttf % num_faces,
                           ttf / num_faces);
      }
    }
  }
  if (set_partition) {
    /* TODO: a copy of this code exists below, make it a function */
    int                 mpirank, mpisize, mpiret;
    int                 first_tree, last_tree, num_trees;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    num_trees = _T8_CMESH_P48_CONN (num_trees);
    first_tree = (mpirank * num_trees) / mpisize;
    last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
#undef _T8_CMESH_P48_CONN
}

t8_cmesh_t
t8_cmesh_new_from_p4est (p4est_connectivity_t * conn, sc_MPI_Comm comm,
                         int do_dup, int set_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 2, comm, do_dup, set_partition);
}

t8_cmesh_t
t8_cmesh_new_from_p8est (p8est_connectivity_t * conn, sc_MPI_Comm comm,
                         int do_dup, int do_partition)
{
  return t8_cmesh_new_from_p4est_ext (conn, 3, comm, do_dup, 1);
}

t8_cmesh_t
t8_cmesh_new_vertex (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[3] = { 0, 0, 0 };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_VERTEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 1);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_line (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[6] = {
    0, 0, 0,
    1, 0, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 2);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_tri (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[9] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 3);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_tet (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    1, 1, 1,
    1, -1, -1,
    -1, 1, -1,
    -1, -1, 1
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TET);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_quad (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[12] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 4);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_hex (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 8);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_pyramid (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[15] = {
    -1, -1, 0,
    1, -1, 0,
    -1, 1, 0,
    1, 1, 0,
    0, 0, sqrt (2)
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PYRAMID);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices,
                              15);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_prism (sc_MPI_Comm comm, int do_dup)
{
  t8_cmesh_t          cmesh;
  double              vertices[18] = {
    0, 0, 0,
    1, 0, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    1, 1, 1
  };

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices, 6);
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_from_class (t8_eclass_t eclass, sc_MPI_Comm comm, int do_dup)
{
  switch (eclass) {
  case T8_ECLASS_VERTEX:
    return t8_cmesh_new_vertex (comm, do_dup);
    break;
  case T8_ECLASS_LINE:
    return t8_cmesh_new_line (comm, do_dup);
    break;
  case T8_ECLASS_TRIANGLE:
    return t8_cmesh_new_tri (comm, do_dup);
    break;
  case T8_ECLASS_QUAD:
    return t8_cmesh_new_quad (comm, do_dup);
    break;
  case T8_ECLASS_TET:
    return t8_cmesh_new_tet (comm, do_dup);
    break;
  case T8_ECLASS_HEX:
    return t8_cmesh_new_hex (comm, do_dup);
    break;
  case T8_ECLASS_PYRAMID:
    return t8_cmesh_new_pyramid (comm, do_dup);
    break;
  case T8_ECLASS_PRISM:
    return t8_cmesh_new_prism (comm, do_dup);
    break;
  default:
    SC_ABORT ("Invalid eclass\n");
    return NULL;
  }
}

/* TODO: This is just a helper function that was needed when we changed the vertex interface
 *       to use attributes. Before we stored a list of vertex coordinates in the cmesh and each tree indexed into this list.
 *       Now each tree carries the coordinates of its vertices.
 *       This function translates from the first approached to the second
 *       and was introduced to avoid rewritting the already existing cmesh_new... functions below.
 *       It would be nice to eventually rewrite these functions correctly.
 */
static void
t8_cmesh_new_translate_vertices_to_attributes (t8_topidx_t * tvertices,
                                               double *vertices,
                                               double *attr_vertices,
                                               int num_vertices)
{
  int                 i;

  for (i = 0; i < num_vertices; i++) {
    attr_vertices[3 * i] = vertices[3 * tvertices[i]];
    attr_vertices[3 * i + 1] = vertices[3 * tvertices[i] + 1];
    attr_vertices[3 * i + 2] = vertices[3 * tvertices[i] + 2];
  }
}

/* The unit cube is constructed from trees of the same eclass.
 * For triangles the square is divided along the (0,0) -- (1,1) diagonal.
 * For prisms the front (y=0) and back (y=1) face are divided into triangles
 * as above.
 */
/* TODO: upgrade with int x,y,z for periodic faces */
t8_cmesh_t
t8_cmesh_new_hypercube (t8_eclass_t eclass, sc_MPI_Comm comm, int do_dup,
                        int do_bcast, int do_partition)
{
  t8_cmesh_t          cmesh;
  int                 num_trees_for_hypercube[T8_ECLASS_COUNT] =
    { 1, 1, 1, 2, 1, 6, 2, 3 };
  int                 i;
  t8_topidx_t         vertices[8];
  double              attr_vertices[24];
  int                 mpirank, mpiret;
  double              vertices_coords[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };

  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  if (!do_bcast || mpirank == 0) {
    t8_cmesh_init (&cmesh);
    for (i = 0; i < num_trees_for_hypercube[eclass]; i++) {
      t8_cmesh_set_tree_class (cmesh, i, eclass);
    }
    switch (eclass) {
    case T8_ECLASS_HEX:
      vertices[4] = 4;
      vertices[5] = 5;
      vertices[6] = 6;
      vertices[7] = 7;
    case T8_ECLASS_QUAD:
      vertices[3] = 3;
      vertices[2] = 2;
    case T8_ECLASS_LINE:
      vertices[1] = 1;
    case T8_ECLASS_VERTEX:
      vertices[0] = 0;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices,
                                                     t8_eclass_num_vertices
                                                     [eclass]);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices,
                                  t8_eclass_num_vertices[eclass]);
      break;
    case T8_ECLASS_PRISM:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 5;
      vertices[3] = 2;
      vertices[4] = 3;
      vertices[5] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      vertices[1] = 5;
      vertices[2] = 4;
      vertices[4] = 7;
      vertices[5] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 6);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 6);
      break;
    case T8_ECLASS_TRIANGLE:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      vertices[0] = 0;
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 3);

      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 3);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 3);
      break;
    case T8_ECLASS_TET:
      t8_cmesh_set_join (cmesh, 0, 1, 1, 2, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 1, 2, 0);
      t8_cmesh_set_join (cmesh, 2, 3, 1, 2, 0);
      t8_cmesh_set_join (cmesh, 3, 4, 1, 2, 0);
      t8_cmesh_set_join (cmesh, 4, 5, 1, 2, 0);
      t8_cmesh_set_join (cmesh, 5, 0, 1, 2, 0);
      vertices[0] = 0;
      vertices[3] = 7;
      vertices[1] = 5;
      vertices[2] = 1;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 1;
      vertices[2] = 3;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 3;
      vertices[2] = 2;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 2;
      vertices[2] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 3, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 6;
      vertices[2] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 4, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      vertices[1] = 4;
      vertices[2] = 5;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 4);
      t8_cmesh_set_tree_vertices (cmesh, 5, t8_get_package_id (), 0,
                                  attr_vertices, 4);
      break;
    case T8_ECLASS_PYRAMID:
      vertices[0] = 1;
      vertices[1] = 3;
      vertices[2] = 0;
      vertices[3] = 2;
      vertices[4] = 7;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      vertices[0] = 0;
      vertices[1] = 2;
      vertices[2] = 4;
      vertices[3] = 6;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 1, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      vertices[0] = 1;
      vertices[1] = 0;
      vertices[2] = 5;
      vertices[3] = 4;
      t8_cmesh_new_translate_vertices_to_attributes (vertices,
                                                     vertices_coords,
                                                     attr_vertices, 5);
      t8_cmesh_set_tree_vertices (cmesh, 2, t8_get_package_id (), 0,
                                  attr_vertices, 5);
      t8_cmesh_set_join (cmesh, 0, 1, 3, 2, 0);
      t8_cmesh_set_join (cmesh, 1, 2, 0, 1, 0);
      t8_cmesh_set_join (cmesh, 2, 0, 2, 0, 0);
      break;
    default:
      break;
    }
  }
  if (do_bcast) {
    if (mpirank != 0) {
      cmesh = NULL;
    }
    cmesh = t8_cmesh_bcast (cmesh, 0, comm);
  }

  if (do_partition) {
    int                 mpirank, mpisize, mpiret;
    int                 first_tree, last_tree, num_trees;

    mpiret = sc_MPI_Comm_rank (comm, &mpirank);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_size (comm, &mpisize);
    SC_CHECK_MPI (mpiret);
    num_trees = num_trees_for_hypercube[eclass];
    first_tree = (mpirank * num_trees) / mpisize;
    last_tree = ((mpirank + 1) * num_trees) / mpisize - 1;
    t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);
  }

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_periodic (sc_MPI_Comm comm, int do_dup, int dim)
{
  t8_cmesh_t          cmesh;
  t8_eclass_t         tree_class;
  double              vertices[24] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1
  };

  T8_ASSERT (dim == 1 || dim == 2 || dim == 3);
  t8_cmesh_init (&cmesh);
  switch (dim) {
  case 1:
    tree_class = T8_ECLASS_LINE;
    break;
  case 2:
    tree_class = T8_ECLASS_QUAD;
    break;
  case 3:
    tree_class = T8_ECLASS_HEX;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  t8_cmesh_set_tree_class (cmesh, 0, tree_class);
  t8_cmesh_set_tree_vertices (cmesh, 0, t8_get_package_id (), 0, vertices,
                              1 << dim);
  t8_cmesh_set_join (cmesh, 0, 0, 0, 1, 0);
  if (dim > 1) {
    t8_cmesh_set_join (cmesh, 0, 0, 2, 3, 0);
  }
  if (dim == 3) {
    t8_cmesh_set_join (cmesh, 0, 0, 4, 5, 0);
  }
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

t8_cmesh_t
t8_cmesh_new_bigmesh (t8_eclass_t eclass, int num_trees, sc_MPI_Comm comm,
                      int do_dup)
{
  t8_cmesh_t          cmesh;
  t8_topidx_t         num_trees_per_eclass[T8_ECLASS_LAST] =
    { 0, 0, 0, 0, 0, 0, 0, 0 };
  int                 i;

  t8_cmesh_init (&cmesh);
  t8_cmesh_set_mpicomm (cmesh, comm, do_dup);
  num_trees_per_eclass[eclass] = num_trees;
  for (i = 0; i < num_trees; i++) {
    t8_cmesh_set_tree_class (cmesh, i, eclass);
  }

  t8_cmesh_commit (cmesh);

  return cmesh;
}
