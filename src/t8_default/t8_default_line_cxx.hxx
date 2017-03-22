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

/** \file t8_default_lines.h
 * The default implementation for lines.
 */

#ifndef T8_DEFAULT_LINE_CXX_HXX
#define T8_DEFAULT_LINE_CXX_HXX

#include <t8_element.h>
#include <t8_element_cxx.hxx>
#include "t8_default_common_cxx.hxx"

T8_EXTERN_C_BEGIN ();

/** Provide an implementation for the line element class.
 * It is written as a self-contained library in the t8_dline_* files.
 */

struct t8_default_scheme_line_c:public t8_default_scheme_common_c
{
public:
  /** The virtual table for a particular implementation of an element class. */

  /** Constructor. */
  t8_default_scheme_line_c ();

  ~t8_default_scheme_line_c ();

/** Return the maximum level allowed for this element class. */
  virtual int         t8_element_maxlevel (void);

/** Return the type of each child in the ordering of the implementation. */
  virtual t8_eclass_t t8_element_child_eclass (int childid)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return T8_ECLASS_ZERO;      /* suppresses compiler warning */
  }

/** Return the refinement level of an element. */
  virtual int         t8_element_level (const t8_element_t * elem);

/** Copy one element to another */
  virtual void        t8_element_copy (const t8_element_t * source,
                                       t8_element_t * dest);

/** Compare to elements. returns negativ if elem1 < elem2, zero if elem1 equals elem2
 *  and positiv if elem1 > elem2.
 *  If elem2 is a copy of elem1 then the elements are equal.
 */
  virtual int         t8_element_compare (const t8_element_t * elem1,
                                          const t8_element_t * elem2)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

/** Construct the parent of a given element. */
  virtual void        t8_element_parent (const t8_element_t * elem,
                                         t8_element_t * parent);

/** Construct a same-size sibling of a given element. */
  virtual void        t8_element_sibling (const t8_element_t * elem,
                                          int sibid, t8_element_t * sibling)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the number of children of an element when it is refined. */
  virtual int         t8_element_num_children (const t8_element_t * elem)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

  /** Return the number of children of an element's face when the element is refined. */
  virtual int         t8_element_num_face_children (const t8_element_t *
                                                    elem, int face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

/** Construct the child element of a given number. */
  virtual void        t8_element_child (const t8_element_t * elem,
                                        int childid, t8_element_t * child);

/** Construct all children of a given element. */
  virtual void        t8_element_children (const t8_element_t * elem,
                                           int length, t8_element_t * c[])
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Return the child id of an element */
  virtual int         t8_element_child_id (const t8_element_t * elem)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

/** Return nonzero if collection of elements is a family */
  virtual int         t8_element_is_family (t8_element_t ** fam)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

/** Construct the nearest common ancestor of two elements in the same tree. */
  virtual void        t8_element_nca (const t8_element_t * elem1,
                                      const t8_element_t * elem2,
                                      t8_element_t * nca)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Return the tree face id given a boundary face. */
  virtual int         t8_element_tree_face (const t8_element_t * elem,
                                            int face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Transform the coordinates of a line considered as boundary element
   *  in a tree-tree connection. */
  virtual void        t8_element_transform_face (const t8_element_t * elem1,
                                                 t8_element_t * elem2,
                                                 int orientation,
                                                 int is_smaller_face);

  /** Given a boundary face inside a root tree's face construct
   *  the element inside the root tree that has the given face as a
   *  face. */
  virtual void        t8_element_extrude_face (const t8_element_t * face,
                                               t8_element_t * elem,
                                               int root_face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Construct the boundary element at a specific face. */
  virtual void        t8_element_boundary_face (const t8_element_t * elem,
                                                int face,
                                                t8_element_t * boundary)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Construct all codimension-one boundary elements of a given element. */
  virtual void        t8_element_boundary (const t8_element_t * elem,
                                           int min_dim, int length,
                                           t8_element_t ** boundary)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

  /** Construct the face neighbor of a given element if this face neighbor
   * is inside the root tree. Return 0 otherwise. */
  virtual int         t8_element_face_neighbor_inside (const t8_element_t *
                                                       elem,
                                                       t8_element_t * neigh,
                                                       int face)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Initialize an element according to a given linear id */
  virtual void        t8_element_set_linear_id (t8_element_t * elem,
                                                int level, uint64_t id);

/** Calculate the linear id of an element */
  virtual u_int64_t   t8_element_get_linear_id (const
                                                t8_element_t *
                                                elem, int level)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

/** Calculate the first descendant of a given element e. That is, the
 *  first element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_first_descendant (const t8_element_t *
                                                   elem, t8_element_t * desc);

/** Calculate the last descendant of a given element e. That is, the
 *  last element in a uniform refinement of e of the maximal possible level.
 */
  virtual void        t8_element_last_descendant (const t8_element_t *
                                                  elem, t8_element_t * desc)
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Compute s as a successor of t*/
  virtual void        t8_element_successor (const t8_element_t * t,
                                            t8_element_t * s, int level);

/** Get the integer coordinates of the anchor node of an element */
  virtual void        t8_element_anchor (const t8_element_t * elem,
                                         int anchor[3])
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

/** Get the integer root length of an element, that is the length of
 *  the level 0 ancestor.
 */
  virtual int         t8_element_root_len (const t8_element_t * elem)
  {
    SC_ABORT ("This function is not implemented yet.\n");
    return 0;                   /* suppresses compiler warning */
  }

  /** Compute the integer coordinates of a given element vertex. */
  virtual void        t8_element_vertex_coords (const t8_element_t * t,
                                                int vertex, int coords[])
  {
    SC_ABORT ("This function is not implemented yet.\n");
  }

};

T8_EXTERN_C_END ();

#endif /* !T8_DEFAULT_LINE_CXX_HXX */