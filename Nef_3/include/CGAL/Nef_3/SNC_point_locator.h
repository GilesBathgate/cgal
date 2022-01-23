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

#ifndef CGAL_NEF_SNC_POINT_LOCATOR_H
#define CGAL_NEF_SNC_POINT_LOCATOR_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_k3_tree_traits.h>
#include <CGAL/Nef_3/K3_tree.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>
#include <string>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 509
#include <CGAL/Nef_2/debug.h>

#undef _CGAL_NEF_TRACEN
#define _CGAL_NEF_TRACEN(msg) CGAL_NEF_TRACEN( "SNC_point_locator: " << msg)

// TODO: find out the proper CGAL replacement for this macro and remove it
#define CGAL_for_each( i, C) for( i = C.begin(); i != C.end(); ++i)

// #define CGAL_NEF_TIMER(instruction) instruction
#define CGAL_NEF_TIMER(instruction)

// #define CGAL_NEF_CLOG(t) std::clog <<" "<<t<<std::endl; std::clog.flush()
#define CGAL_NEF_CLOG(t)

namespace CGAL {

template <typename SNC_decorator>
class SNC_point_locator
{
 public:
  class Intersection_call_back;
  typedef SNC_decorator Base;
  typedef SNC_point_locator<SNC_decorator> Self;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
protected:
  std::string version_;
  // time for construction, point location, ray shooting and intersection test
  mutable Timer ct_t, pl_t, rs_t, it_t;

public:
  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Aff_transformation_3
                                  Aff_transformation_3;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::Volume_handle Volume_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;


  const std::string& version() const { return version_; }

  virtual Object_handle locate(const Point_3& p) const = 0;

  virtual Object_handle shoot(const Ray_3& s, int mask=255) const = 0;

  virtual void intersect_with_edges( Halfedge_handle edge,
                                     const Intersection_call_back& call_back)
    const = 0;

  virtual void intersect_with_facets( Halfedge_handle edge,
                                      const Intersection_call_back& call_back)
    const = 0;

  virtual void intersect_with_edges_and_facets( Halfedge_handle edge,
        const Intersection_call_back& call_back) const = 0;

  class Intersection_call_back
  {
  public:
    virtual void operator()( Halfedge_handle edge, Object_handle object,
                             const Point_3& intersection_point) const = 0;

    virtual ~Intersection_call_back() {}
  };

  virtual void initialize(SNC_structure* W) = 0;

  virtual Self* clone() const = 0;

  virtual void transform(const Aff_transformation_3& t) = 0;

  virtual void add_facet(Halffacet_handle) {}

  virtual void add_edge(Halfedge_handle) {}

  virtual void add_vertex(Vertex_handle) {}

  virtual ~SNC_point_locator() noexcept(!CGAL_ASSERTIONS_ENABLED)
  {
    CGAL_NEF_CLOG("");
    CGAL_NEF_CLOG("construction_time:  "<<ct_t.time());
    CGAL_NEF_CLOG("pointlocation_time: "<<pl_t.time());
    CGAL_NEF_CLOG("rayshooting_time:   "<<rs_t.time());
    CGAL_NEF_CLOG("intersection_time:  "<<it_t.time());
    // warning: the total time showed here could be actually larger
    // that the real time used by this class, since point location
    // and intersection test use the ray shooter and so the same time
    // could be account to more than one timer
    CGAL_NEF_CLOG("pltotal_time:       "<<
      ct_t.time()+pl_t.time()+rs_t.time()+it_t.time());
  };
};

template <typename SNC_decorator>
class SNC_point_locator_by_spatial_subdivision :
  public SNC_point_locator<SNC_decorator>,
  public SNC_decorator
{
  typedef SNC_decorator Base;
  typedef CGAL::SNC_point_locator<SNC_decorator> SNC_point_locator;
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_decorator> Self;
  typedef typename SNC_decorator::SNC_structure SNC_structure;
  typedef typename SNC_decorator::Decorator_traits Decorator_traits;
  typedef typename Decorator_traits::SM_decorator SM_decorator;
  typedef CGAL::SNC_intersection<SNC_structure> SNC_intersection;
  typedef typename SNC_decorator::SM_point_locator SM_point_locator;
  typedef typename SNC_decorator::Kernel Kernel;

public:
  typedef typename CGAL::SNC_k3_tree_traits<SNC_decorator> K3_tree_traits;

  typedef typename CGAL::K3_tree<K3_tree_traits> K3_tree;
  typedef typename K3_tree::Node_handle Node_handle;
  typedef K3_tree SNC_candidate_provider;

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Triangle_3 Triangle_3;
  typedef typename SNC_structure::Aff_transformation_3
                                  Aff_transformation_3;

  typedef typename SNC_structure::Infi_box Infi_box;

  typedef typename Decorator_traits::Vertex_handle Vertex_handle;
  typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
  typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
  typedef typename Decorator_traits::SHalfedge_handle SHalfedge_handle;
  typedef typename Decorator_traits::SFace_handle SFace_handle;
  typedef typename Decorator_traits::Vertex_iterator Vertex_iterator;
  typedef typename Decorator_traits::Halfedge_iterator Halfedge_iterator;
  typedef typename Decorator_traits::Halffacet_iterator Halffacet_iterator;
  typedef typename Decorator_traits::Volume_handle Volume_handle;

  typedef typename Decorator_traits::Halffacet_cycle_iterator
                                     Halffacet_cycle_iterator;
  typedef typename Decorator_traits::SHalfedge_around_facet_circulator
                                     SHalfedge_around_facet_circulator;

  typedef typename SNC_candidate_provider::Object_list Object_list;
  typedef typename Object_list::iterator Object_list_iterator;
  typedef typename SNC_candidate_provider::Vertex_const_iterator Vertex_list_iterator;
  typedef typename SNC_candidate_provider::Halfedge_const_iterator Halfedge_list_iterator;
  typedef typename SNC_candidate_provider::Halffacet_const_iterator Halffacet_list_iterator;

  typedef typename SNC_candidate_provider::Objects_along_ray Objects_along_ray;
  typedef typename Objects_along_ray::Iterator Objects_along_ray_iterator;

  using Base::get_visible_facet;
public:
  SNC_point_locator_by_spatial_subdivision() :
    initialized(false), candidate_provider(0) {}


  virtual void initialize(SNC_structure* W) {

    if(initialized)
      delete candidate_provider;

    this->set_snc(*W);
    candidate_provider = new SNC_candidate_provider(W);

    initialized = true;
  }

  virtual Self* clone() const {
    return new Self;
  }

  virtual void transform(const Aff_transformation_3& t) {
    candidate_provider->transform(t);
  }

  virtual ~SNC_point_locator_by_spatial_subdivision() noexcept(!CGAL_ASSERTIONS_ENABLED)
  {
    CGAL_destructor_warning(initialized ||
                 candidate_provider == 0); // required?
    if(initialized)
      delete candidate_provider;
  }

  static inline bool v_vertex_of_e(Vertex_handle v, Halfedge_handle e) {
    return (e->source() == v) || (e->twin()->source() == v);
  }

  // We next check if v is a vertex on the face to avoid a geometric test
  static inline bool v_vertex_of_f(Vertex_handle v, Halffacet_handle f) {
    Halffacet_cycle_iterator fci;
    for(fci=f->facet_cycles_begin(); fci!=f->facet_cycles_end(); ++fci) {
      if(fci.is_shalfedge()) {
        SHalfedge_around_facet_circulator sfc(fci), send(sfc);
        CGAL_For_all(sfc,send) {
          if(sfc->source()->center_vertex() ==  v){
            return true;
          }
        }
      }
    }
    return false;
  }

  virtual Object_handle shoot(const Ray_3& ray, int mask=255) const {
    CGAL_NEF_TIMER(rs_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "shooting: "<<ray);

    Vertex_handle v_result;
    Halfedge_handle e_result;
    Halffacet_handle f_result;
    enum What { None=false, Vertex, Edge, Facet };
    What hit = None;
    Point_3 eor = CGAL::ORIGIN; // 'end of ray', the latest ray's hit point
    Objects_along_ray objects = candidate_provider->objects_along_ray(ray);
    for(Objects_along_ray_iterator objects_iterator = objects.begin();
        objects_iterator != objects.end(); ++objects_iterator) {
      Node_handle node = objects_iterator.get_node();
      if(((mask&1) != 0)) {
        for(Vertex_list_iterator vi=node->vertices_begin(); vi!=node->vertices_end(); ++vi) {
          Vertex_handle v = *vi;
          _CGAL_NEF_TRACEN("trying vertex on "<<v->point());
          if( ray.source() != v->point() && ray.has_on(v->point())) {
            _CGAL_NEF_TRACEN("the ray intersects the vertex");
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            CGAL_assertion_code
              (if( hit )_CGAL_NEF_TRACEN("prev. intersection on "<<eor));
            if( hit && !Segment_3( ray.source(), eor).has_on(v->point()))
              continue;
            eor = v->point();
            v_result = v;
            hit = Vertex;
            _CGAL_NEF_TRACEN("the vertex becomes the new hit object");
          }
        }
      }
      if(((mask&2) != 0)) {
        for(Halfedge_list_iterator ei=node->edges_begin(); ei!=node->edges_end(); ++ei) {
          Halfedge_handle e = *ei;
          Point_3 q;
          _CGAL_NEF_TRACEN("trying edge on "<< Segment_3(e->source()->point(),e->twin()->source()->point()));
          if( !v_vertex_of_e(v_result,e) && is.does_intersect_internally( ray, Segment_3(e->source()->point(),
                                                           e->twin()->source()->point()), q)) {
            _CGAL_NEF_TRACEN("ray intersects edge on "<<q);
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            CGAL_assertion_code
              (if( hit) _CGAL_NEF_TRACEN("prev. intersection on "<<eor));
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _CGAL_NEF_TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
                continue;
            eor = q;
            e_result = e;
            hit = Edge;
            _CGAL_NEF_TRACEN("the edge becomes the new hit object");
        }
      }
      }
      if(((mask&4) != 0)) {
        for(Halffacet_list_iterator fi=node->facets_begin(); fi!=node->facets_end(); ++fi) {
          Halffacet_handle f = *fi;
          Point_3 q;
          _CGAL_NEF_TRACEN("trying facet with on plane "<<f->plane()<<
                  " with point on "<<f->plane().point());
          if( !v_vertex_of_f(v_result,f) && is.does_intersect_internally( ray, f, q, true) ) {
            _CGAL_NEF_TRACEN("ray intersects facet on "<<q);
            _CGAL_NEF_TRACEN("prev. intersection? "<<hit);
            if( hit) { _CGAL_NEF_TRACEN("prev. intersection on "<<eor); }
            if( hit && !has_smaller_distance_to_point( ray.source(), q, eor))
              continue;
            _CGAL_NEF_TRACEN("is the intersection point on the current cell? "<<
                    candidate_provider->is_point_on_cell( q, objects_iterator));
            if( !candidate_provider->is_point_on_cell( q, objects_iterator))
              continue;
            eor = q;
            f_result = f;
            hit = Facet;
            _CGAL_NEF_TRACEN("the facet becomes the new hit object");
          }
        }
      }
      if(hit) break;
    }
    CGAL_NEF_TIMER(rs_t.stop());
    switch(hit) {
      case Facet: return make_object(f_result);
      case Edge: return make_object(e_result);
      case Vertex: return make_object(v_result);
      default: return Object_handle();
    }
  }

  virtual Object_handle locate( const Point_3& p) const {
    if(Infi_box::extended_kernel()) {
    CGAL_NEF_TIMER(pl_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "locate "<<p);

    Node_handle candidates = candidate_provider->nodes_around_point(p);

    for(Vertex_list_iterator vi = candidates->vertices_begin();
        vi != candidates->vertices_end(); ++vi) {
      Vertex_handle v = *vi;
      if ( p == v->point()) {
        _CGAL_NEF_TRACEN("found on vertex "<<v->point());
        CGAL_NEF_TIMER(pl_t.stop());
        return make_object(v);
      }
    }
    for(Halfedge_list_iterator ei = candidates->edges_begin();
        ei != candidates->edges_end(); ++ei) {
      Halfedge_handle e = *ei;
      if ( is.does_contain_internally(Segment_3(e->source()->point(),e->twin()->source()->point()), p) ) {
        _CGAL_NEF_TRACEN("found on edge "<<Segment_3(e->source()->point(),e->twin()->source()->point()));
        CGAL_NEF_TIMER(pl_t.stop());
        return make_object(e);
      }
    }
    for(Halffacet_list_iterator fi = candidates->facets_begin();
        fi != candidates->facets_end(); ++fi) {
      Halffacet_handle f = *fi;
      if (is.does_contain_internally( f, p) ) {
        _CGAL_NEF_TRACEN("found on facet...");
        CGAL_NEF_TIMER(pl_t.stop());
        return make_object(f);
      }
    }

    _CGAL_NEF_TRACEN("point not found in 2-skeleton");
    _CGAL_NEF_TRACEN("shooting ray to determine the volume");
    Ray_3 r( p, Vector_3( -1, 0, 0));
    CGAL_NEF_TIMER(pl_t.stop());
    return make_object(determine_volume(r));

  } else {   // standard kernel


    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "locate "<<p);
    typedef typename SNC_structure::FT FT;

    Node_handle candidates = candidate_provider->nodes_around_point(p);

    if(candidates->empty())
      return make_object(Base(*this).volumes_begin());

    Vertex_list_iterator vi = candidates->vertices_begin();

    Vertex_handle v_result = *vi;
    if(p==v_result->point())
      return make_object(v_result);

    FT min_distance = CGAL::squared_distance(v_result->point(),p);
    ++vi;

    for(;vi!=candidates->vertices_end(); ++vi) {
      Vertex_handle v = *vi;
      if ( p == v->point()) {
        _CGAL_NEF_TRACEN("found on vertex "<<v->point());
        return make_object(v);
      }
      FT tmp_distance = CGAL::squared_distance(v->point(),p);
      if(tmp_distance < min_distance) {
        min_distance = tmp_distance;
        v_result = v;
      }
    }

    Segment_3 s(p,v_result->point());

    Halfedge_handle e_result;
    for(Halfedge_list_iterator ei=candidates->edges_begin(); ei!=candidates->edges_end(); ++ei) {
      Halfedge_handle e = *ei;
      Segment_3 ss(e->source()->point(),e->twin()->source()->point());
      CGAL_NEF_TRACEN("test edge " << e->source()->point() << "->" << e->twin()->source()->point());
      if (is.does_contain_internally(ss, p)) {
      _CGAL_NEF_TRACEN("found on edge "<< ss);
        return make_object(e);
      }
      Point_3 ip;
      if(!v_vertex_of_e(v_result,e) && is.does_intersect_internally(s, ss, ip)) {
        s = Segment_3(p, normalized(ip));
        e_result = e;
      }
    }

    Halffacet_handle f_result;
    for(Halffacet_list_iterator fi=candidates->facets_begin(); fi!=candidates->facets_end(); ++fi) {
      Halffacet_handle f = *fi;
      CGAL_NEF_TRACEN("test facet " << f->plane());
      if (is.does_contain_internally(f,p) ) {
        _CGAL_NEF_TRACEN("found on facet...");
        return make_object(f);
      }
      Point_3 ip;
      if( !v_vertex_of_f(v_result,f) &&  is.does_intersect_internally(s,f,ip) ) {
        s = Segment_3(p, normalized(ip));
        f_result = f;
      }
    }

    if( f_result != nullptr ) {
      _CGAL_NEF_TRACEN("facet hit, obtaining volume...");
      if(f_result->plane().oriented_side(p) == ON_NEGATIVE_SIDE)
        f_result = f_result->twin();
      return make_object(f_result->incident_volume());
    } else if( e_result != nullptr ) {
      SM_decorator SD(&*e_result->source());
      if( SD.is_isolated(e_result))
        return make_object(e_result->incident_sface()->volume());
      return make_object(get_visible_facet(e_result,Ray_3(s.source(),s.to_vector()))->incident_volume());
    }

    _CGAL_NEF_TRACEN("vertex hit, obtaining volume..." << v_result->point());
    SM_point_locator L(&*v_result);
    Object_handle so = L.locate(s.source()-s.target());
    SFace_handle sf;
    if(CGAL::assign(sf,so))
      return make_object(sf->volume());

    CGAL_error_msg( "wrong handle type");
    return Object_handle();
  }
  }

  virtual void intersect_with_edges_and_facets( Halfedge_handle e0,
        const typename SNC_point_locator::Intersection_call_back& call_back) const {

    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));


    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, Segment_3(e->source()->point(),
                                                       e->twin()->source()->point()), q)) {
          q = normalized(q);
          call_back( e0, make_object(Halfedge_handle(e)), q);
          _CGAL_NEF_TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
                                                                e->twin()->source()->point())<<" on "<<q);
        }
      }
      else if( CGAL::assign( f, *o)) {
#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, f, q) ) {
          q = normalized(q);
          call_back( e0, make_object(Halffacet_handle(f)), q);
          _CGAL_NEF_TRACEN("edge intersects facet on plane "<<f->plane()<<" on "<<q);
        }
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

  virtual void intersect_with_edges( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<<&*e0<<' '<<Segment_3(e0->source()->point(),
                                                         e0->twin()->source()->point()));
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, Segment_3(e->source()->point(),
                                                       e->twin()->source()->point()), q)) {
          q = normalized(q);
          call_back( e0, make_object(Halfedge_handle(e)), q);
          _CGAL_NEF_TRACEN("edge intersects edge "<<' '<<&*e<< Segment_3(e->source()->point(),
                                                                e->twin()->source()->point())<<" on "<<q);
        }
      }
      else if( CGAL::assign( f, *o)) {
        /* do nothing */
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

  virtual void intersect_with_facets( Halfedge_handle e0,
    const typename SNC_point_locator::Intersection_call_back& call_back) const {
    CGAL_NEF_TIMER(it_t.start());
    CGAL_assertion( initialized);
    _CGAL_NEF_TRACEN( "intersecting edge: "<< Segment_3(e0->source()->point(),
                                               e0->twin()->source()->point()));
    Segment_3 s(Segment_3(e0->source()->point(),e0->twin()->source()->point()));
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Object_list_iterator o;
    Object_list objects = candidate_provider->objects_around_segment(s);
    CGAL_for_each( o, objects) {
      if( CGAL::assign( v, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( e, *o)) {
        /* do nothing */
      }
      else if( CGAL::assign( f, *o)) {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersection_candidates;
#endif

        Point_3 q;
        if( is.does_intersect_internally( s, f, q) ) {
          q = normalized(q);
          call_back( e0, make_object(Halffacet_handle(f)), q);
          _CGAL_NEF_TRACEN("edge intersects facet on plane "<<f->plane()<<" on "<<q);
        }
      }
      else
        CGAL_error_msg( "wrong handle");
    }
    CGAL_NEF_TIMER(it_t.stop());
  }

private:
  Volume_handle determine_volume( const Ray_3& ray) const {
    Halffacet_handle f_below;
    Object_handle o = shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    if( CGAL::assign( v, o)) {
      _CGAL_NEF_TRACEN("vertex hit, obtaining volume...");
      f_below = get_visible_facet( v, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      SM_decorator SD(&*v); // now, the vertex has no incident facets
      CGAL_assertion( SD.number_of_sfaces() == 1);
      return SD.sfaces_begin()->volume();
    }
    else if( CGAL::assign( e, o)) {
      _CGAL_NEF_TRACEN("edge hit, obtaining volume...");
      f_below = get_visible_facet( e, ray);
      if( f_below != Halffacet_handle())
        return f_below->incident_volume();
      CGAL_assertion_code(SM_decorator SD(&*e->source())); // now, the edge has no incident facets
      CGAL_assertion(SD.is_isolated(e));
      return e->incident_sface()->volume();
    }
    else if( CGAL::assign( f, o)) {
      _CGAL_NEF_TRACEN("facet hit, obtaining volume...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
      return f_below->incident_volume();
    }
    return Base(*this).volumes_begin(); // TODO: Comment this hack!
  }

public:
  void add_facet(Halffacet_handle f) {
    candidate_provider->add_facet(f);
  }

  void add_edge(Halfedge_handle e) {
    candidate_provider->add_edge(e);
  }

  void add_vertex(Vertex_handle v) {
    candidate_provider->add_vertex(v);
  }

private:
  bool initialized;
  SNC_candidate_provider* candidate_provider;
  SNC_intersection is;
};


} //namespace CGAL
#endif // CGAL_NEF_SNC_POINT_LOCATOR_H
