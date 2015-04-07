/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2010 The University of Texas System
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

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

#include "t8_default_common.h"
#include "t8_default_tet.h"

typedef int8_t      t8_default_tet_type_t;
typedef int8_t      t8_default_tet_cube_id_t;

struct t8_default_tet_id
{
  t8_default_tet_type_t type;
  t8_tcoord_t         anchor_coordinates[3];
};

int                 t8_tet_cid_type_to_parenttype[8][6] = {
  {0, 1, 2, 3, 4, 5},
  {0, 1, 1, 1, 0, 0},
  {2, 2, 2, 3, 3, 3},
  {1, 1, 2, 2, 2, 1},
  {5, 5, 4, 4, 4, 5},
  {0, 0, 0, 5, 5, 5},
  {4, 3, 3, 3, 4, 4},
  {0, 1, 2, 3, 4, 5}
};

/* In dependence of a type x give the type of
 * the child with Bey number y */
int                 t8_tet_type_of_child[6][8] = {
  {0, 0, 0, 0, 4, 5, 2, 1},
  {1, 1, 1, 1, 3, 2, 5, 0},
  {2, 2, 2, 2, 0, 1, 4, 3},
  {3, 3, 3, 3, 5, 4, 1, 2},
  {4, 4, 4, 4, 2, 3, 0, 5},
  {5, 5, 5, 5, 1, 0, 3, 4}
};

static              size_t
t8_default_tet_size (void)
{
  return sizeof (t8_tet_t);
}

static t8_default_tet_id_t *
t8_default_tet_id_new (t8_default_tet_type_t type,
                       t8_tcoord_t anchor_coordinates[3])
{
  t8_default_tet_id_t *tet_id;
  int                 i;

  tet_id = T8_ALLOC (t8_default_tet_id_t, 1);
  tet_id->type = type;
  for (i = 0; i < 3; i++) {
    tet_id->anchor_coordinates[i] = anchor_coordinates[i];
  }
  return tet_id;
}

static void
t8_default_tet_id_destroy (t8_default_tet_id_t * tid)
{
  T8_FREE (tid);
}

static              t8_default_tet_type_t
t8_default_tet_get_type (const t8_tet_t * t)
{
  return t->tet_id->type;
}

static void
t8_default_tet_set_type (t8_tet_t * t, t8_default_tet_type_t type)
{
  T8_ASSERT (0 <= type && type <= 5);

  t->tet_id->type = type;
}

static              t8_tcoord_t
t8_default_tet_get_coordinate (const t8_tet_t * t, int i)
{
  T8_ASSERT (0 <= i && i < 3);
  return t->tet_id->anchor_coordinates[i];
}

static void
t8_default_tet_set_coordinate (t8_tet_t * t, int i, t8_tcoord_t value)
{
  T8_ASSERT (0 <= i && i < 3);
  t->tet_id->anchor_coordinates[i] = value;
}

static              t8_default_tet_cube_id_t
t8_default_tet_compute_cubeid (const t8_default_tet_id_t * tid, int level)
{
  t8_default_tet_cube_id_t id = 0;
  t8_tcoord_t         h;

  T8_ASSERT (0 <= level && level <= T8_TET_MAX_LEVEL);
  h = T8_TET_ROOT_LEN (level);

  if (level == 0) {
    return 0;
  }

  id |= ((tid->anchor_coordinates[0] & h) ? 0x01 : 0);
  id |= ((tid->anchor_coordinates[1] & h) ? 0x02 : 0);
  id |= ((tid->anchor_coordinates[2] & h) ? 0x04 : 0);

  return id;
}

static t8_default_tet_id_t *
t8_default_tet_parent_tetid (const t8_default_tet_id_t * tid,
                             const int8_t level)
{
  t8_default_tet_cube_id_t cid;
  t8_default_tet_type_t parent_type;
  t8_tcoord_t         parent_coord[3], h;
  int                 i;

  /* Compute type of parent */
  cid = t8_default_tet_compute_cubeid (tid, level);
  parent_type = t8_tet_cid_type_to_parenttype[cid][tid->type];
  /* Compute anchor coordinates of parent */
  h = T8_TET_ROOT_LEN (level);
  for (i = 0; i < 3; i++) {
    parent_coord[i] = tid->anchor_coordinates[i] & ~h;
  }

  return t8_default_tet_id_new (parent_type, parent_coord);
}

static void
t8_default_tet_parent (const t8_element_t * elem, t8_element_t * parent)
{
  const t8_tet_t     *t = (const t8_tet_t *) elem;
  t8_tet_t           *p = (t8_tet_t *) parent;
  t8_default_tet_cube_id_t cid;
  t8_tcoord_t         h;
  int                 i;

  T8_ASSERT (t->level > 0);

  p->itype = t->itype;
  p->level = t->level - 1;

  /* Compute type of parent */
  cid = t8_default_tet_compute_cubeid (t->tet_id, t->level);
  t8_default_tet_set_type (p, t8_tet_cid_type_to_parenttype[cid]
                           [t8_default_tet_get_type (t)]);
  /* Set coordinates of parent */
  h = T8_TET_ROOT_LEN (t->level);
  for (i = 0; i < 3; i++) {
    p->tet_id->anchor_coordinates[i] = t->tet_id->anchor_coordinates[i] & ~h;
  }
}

static void
t8_default_tet_compute_coords (const t8_tet_t * t,
                               t8_tcoord_t coordinates[4][3])
{
  t8_default_tet_type_t type;
  int                 ei, ej;
  int                 i, j;
  t8_tcoord_t         h;

  type = t8_default_tet_get_type (t);
  h = T8_TET_ROOT_LEN (t->level);
  ei = type / 2;
  if (type % 2 == 0) {
    ej = (ei + 2) % 3;
  }
  else {
    ej = (ei + 1) % 3;
  }
  for (i = 0; i < 2; i++) {
    coordinates[0][i] = t8_default_tet_get_coordinate (t, i);
    coordinates[1][i] = coordinates[0][i];
    coordinates[2][i] = coordinates[0][i];
    coordinates[3][i] = coordinates[0][i] + h;
  }
  coordinates[1][ei] += h;
  coordinates[2][ei] += h;
  coordinates[2][ej] += h;
}

/* The childid here is the Bey child id,
 * not the Morton child id */
static void
t8_default_tet_child (const t8_element_t * elem, int childid,
                      t8_element_t * child)
{
  const t8_tet_t     *t = (const t8_tet_t *) (elem);
  t8_tet_t           *c = (t8_tet_t *) (child);
  t8_tcoord_t         t_coordinates[4][3], temp_coord;
  t8_default_tet_type_t type;
  int                 coord2, i;

  T8_ASSERT (t->level < T8_TET_MAX_LEVEL);
  T8_ASSERT (0 <= childid && childid < 8);

  /* Compute anchor coordinates of child */
  if (childid == 0) {
    for (i = 0; i < 3; i++) {
      t8_default_tet_set_coordinate (c, i,
                                     t8_default_tet_get_coordinate (t, i));
    }
  }
  else {
    switch (childid) {
    case 1:
    case 4:
    case 5:
      coord2 = 1;
      break;
    case 2:
    case 6:
    case 7:
      coord2 = 2;
      break;
    case 3:
      coord2 = 3;
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
    /* i-th anchor coordinate of child is (X_(0,i)+X_(coord2,i))/2
     * where X_(i,j) is the j-th coordinate of t's ith node */
    t8_default_tet_compute_coords (t, t_coordinates);
    for (i = 0; i < 3; i++) {
      temp_coord = (t_coordinates[0][i] + t_coordinates[coord2][i]) >> 1;
      t8_default_tet_set_coordinate (c, i, temp_coord);
    }
  }

  /* Compute type of child */
  type = t8_default_tet_get_type (t);
  t8_default_tet_set_type (c, t8_tet_type_of_child[type][childid]);

  c->level = t->level + 1;
}

t8_type_scheme_t   *
t8_default_scheme_new_tet (void)
{
  t8_type_scheme_t   *ts;

  ts = T8_ALLOC (t8_type_scheme_t, 1);

  ts->elem_size = t8_default_tet_size;
  ts->elem_parent = t8_default_tet_parent;
  ts->elem_child = t8_default_tet_child;

  ts->elem_new = t8_default_mempool_alloc;
  ts->elem_destroy = t8_default_mempool_free;
  ts->ts_destroy = t8_default_scheme_mempool_destroy;
  ts->ts_context = sc_mempool_new (sizeof (t8_tet_t));

  return ts;
}