// Copyright (c) 2020  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_PMP_IO_POLYGON_MESH_IO_H
#define CGAL_PMP_IO_POLYGON_MESH_IO_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
  \ingroup PMP_IO_grp

 * \brief Attempts to read a file as a polygon mesh; in case of failure, reads the file as a polygon soup,
 * repairs and orients it to obtain a polygon mesh.
 *
 * Supported file formats are the following:
 * - \ref IOStreamOFF (`.off`)
 * - \ref IOStreamOBJ (`.obj`)
 * - \ref IOStreamSTL (`.stl`)
 * - \ref IOStreamPLY (`.ply`)
 * - \ref IOStreamGocad (`.ts`)
 * - \ref IOStreamVTK (`.vtp`)
 *
 * The format is detected from the filename extension.
 *
 * \tparam PolygonMesh a model of `MutableFaceGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param fname the name of the input file
 * \param g the polygon mesh
 * \param np sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `g`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, g)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{repair_polygon_soup}
 *     \cgalParamDescription{a parameter used indicate whether `CGAL::Polygon_mesh_processing::repair_polygon_soup()`
 *                           should be called on the soup in case of issues in the input.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{whether extra information is printed when an incident occurs during reading}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`true`}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \return `true` if the reading and conversion worked, `false` otherwise.
 *
 * \sa \link PkgBGLIOFct `CGAL::write_polygon_mesh()` \endlink
 */
template <typename PolygonMesh, typename NamedParameters>
bool read_polygon_mesh(const char* fname,
                       PolygonMesh& g,
                       const NamedParameters& np)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type                     Point;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool ok = CGAL::read_polygon_mesh(fname, g, np);

  if(ok)
    return true;

  clear(g);

  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  if(!CGAL::read_polygon_soup(fname, points, faces))
  {
    std::cerr << "Error: cannot read file\n";
    return false;
  }

  std::cout << "Cleaning polygon soup..." << std::endl;
  const bool do_repair = choose_parameter(get_parameter(np, internal_np::repair_polygon_soup), true);
  if(do_repair)
    PMP::repair_polygon_soup(points, faces, np);

  if(!PMP::orient_polygon_soup(points, faces))
    std::cerr << "W: File does not describe a polygon mesh" << std::endl;

  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    return false;

  PMP::polygon_soup_to_polygon_mesh(points, faces, g, parameters::all_default(), np);

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename PolygonMesh>
bool read_polygon_mesh(const char* fname, PolygonMesh& g)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname, g, parameters::all_default());
}

template <typename PolygonMesh, typename NamedParameters>
bool read_polygon_mesh(const std::string& fname, PolygonMesh& g, const NamedParameters& np)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname.c_str(), g, np);
}

template <typename PolygonMesh>
bool read_polygon_mesh(const std::string& fname, PolygonMesh& g)
{
  return CGAL::Polygon_mesh_processing::read_polygon_mesh(fname, g, parameters::all_default());
}

/// \endcond

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_PMP_IO_POLYGON_MESH_IO_H