// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SNC_K3_TREE_TRAITS_H
#define CGAL_NEF_SNC_K3_TREE_TRAITS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_3/Bounding_box_3.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Unique_hash_map.h>
#include <list>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 503
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class SNC_decorator>
class Side_of_plane {

  typedef typename SNC_decorator::Decorator_traits Decorator_traits;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;

  typedef typename Decorator_traits::Halffacet_cycle_iterator
    Halffacet_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator
    SHalfedge_around_facet_circulator;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;

  typedef typename SNC_decorator::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;

  static inline constexpr Oriented_side unknown_side = static_cast<Oriented_side>(-2);

public:
  Side_of_plane(const Point_3& p, int c) : OnSideMap(unknown_side), coord(c), pop(p) {
    CGAL_assertion( c >= 0 && c <= 2 );
  }
  void reserve(std::size_t n) { OnSideMap.reserve(n); }

  Oriented_side operator()(Vertex_handle v) {
  Comparison_result cr;
    Oriented_side& side = OnSideMap[v];
    if(side == unknown_side) {
      cr = (coord == 0) ? CGAL::compare_x(v->point(), pop) :
           (coord == 1) ? CGAL::compare_y(v->point(), pop) :
         /*(coord == 2)*/ CGAL::compare_z(v->point(), pop);

      side = cr == LARGER ? ON_POSITIVE_SIDE :
             cr == SMALLER ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
    }
    return side;
  }

  /*
   An edge is considered intersecting a plane if its endpoints lie on the
   plane or if they lie on different sides.  Partial tangency is not considered
   as intersection, due the fact that a lower dimensional face (the vertex)
   should be already reported as an object intersecting the plane.
  */
  Oriented_side operator()(Halfedge_handle e) {
    Vertex_handle vs = e->source();
    Vertex_handle vt = e->twin()->source();

    Oriented_side src_side = (*this) (vs);
    Oriented_side tgt_side = (*this) (vt);
    if( src_side == tgt_side)
      return src_side;
    if( src_side == ON_ORIENTED_BOUNDARY)
      return tgt_side;
    if( tgt_side == ON_ORIENTED_BOUNDARY)
      return src_side;
    return ON_ORIENTED_BOUNDARY;
  }

  /*
   As for the edges, if a facet is tangent to the plane it is not considered as
   a intersection since lower dimensional faces, like the edges and vertices
   where the tangency occurs, should be reported as the objects intersecting
   the plane.
   A valid flat facet requires only three vertices to establish coplanarity.
   However, to safely handle arbitrary or degenerate topologies, the loop
   completely traverses the cycle so the intersection is known as far as two
   vertices are located on different sides of the plane.
  */
  Oriented_side operator()(Halffacet_handle f) {
    CGAL_assertion(f->facet_cycles_begin() != f->facet_cycles_end());

    Halffacet_cycle_iterator fc(f->facet_cycles_begin());
    CGAL_assertion(fc.is_shalfedge());
    SHalfedge_handle e(fc);
    SHalfedge_around_facet_circulator sc(e), send(sc);
    //CGAL_assertion( iterator_distance( sc, send) >= 3); // TODO: facet with 2 vertices was found, is it possible?

    bool all_on_boundary = true;
    Oriented_side facet_side;
    do {
      const Vertex_handle& v = sc->source()->center_vertex();
      facet_side = (*this)(v);
      if (facet_side != ON_ORIENTED_BOUNDARY) {
        all_on_boundary = false;
        break;
      }
      ++sc;
    } while (sc != send);

    if (all_on_boundary)
      return ON_ORIENTED_BOUNDARY;
    CGAL_assertion( facet_side != ON_ORIENTED_BOUNDARY);
    Oriented_side point_side;
    ++sc;
    while (sc != send) {
      const Vertex_handle& v = sc->source()->center_vertex();
      point_side = (*this)(v);
      if (point_side != ON_ORIENTED_BOUNDARY && point_side != facet_side)
        return ON_ORIENTED_BOUNDARY;
      ++sc;
    }

    return facet_side;
  }
private:
  Unique_hash_map<Vertex_handle,Oriented_side> OnSideMap;
  int coord;
  const Point_3 pop;
};

template <class Decorator>
class SNC_k3_tree_traits {

public:
  typedef Decorator SNC_decorator;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
  typedef typename SNC_structure::Kernel Kernel;

  typedef typename SNC_structure::Infi_box Infimaximal_box;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef std::vector<Object_handle> Object_list;
  typedef std::vector<Vertex_handle> Vertex_list;
  typedef std::vector<Halfedge_handle> Halfedge_list;
  typedef std::vector<Halffacet_handle> Halffacet_list;

  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Ray_3 Ray_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;

  typedef typename Kernel::RT RT;
  typedef typename Kernel::Kernel_tag Kernel_tag;
  typedef CGAL::Bounding_box_3<Tag_true, Kernel> Bounding_box_3;
  typedef typename Kernel::Intersect_3 Intersect;
  typedef CGAL::Side_of_plane<SNC_decorator> Side_of_plane;

  Intersect intersect_object() const {
    return Intersect();
  }

};

} //namespace CGAL

#endif // CGAL_NEF_SNC_K3_TREE_TRAITS_H
