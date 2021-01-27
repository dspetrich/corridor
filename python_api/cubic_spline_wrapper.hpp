#pragma once

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numpy.hpp>
#include <utility>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"
#include "utility.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace cr = corridor;
namespace cs = cr::cubic_spline;

// /////////////////////////////////////////////////////////////////////////////
// Cubic Spline wrapper
// /////////////////////////////////////////////////////////////////////////////

struct CubicSplineWrapper {
  // CPP object
  cs::CubicSpline cubic_spline_;

  CubicSplineWrapper(const int id, const py::list& node_x,
                     const py::list& node_y) {
    using namespace corridor;
    CartesianPoints2D points = toCartesianPoints(node_x, node_y);
    cubic_spline_ = cs::CubicSpline(id, points);
  }

  CubicSplineWrapper(const int id, const py::list& node_x,
                     const py::list& node_y, const py::list& py_first_tangent,
                     const py::list& py_last_tangent) {
    using namespace corridor;
    CartesianPoints2D points = toCartesianPoints(node_x, node_y);
    CartesianVector2D first_tangent = to_cartesian_vector(py_first_tangent);
    CartesianVector2D last_tangent = to_cartesian_vector(py_last_tangent);
    cubic_spline_ = cs::CubicSpline(id, points, first_tangent, last_tangent);
  }
  const cr::RealType GetTotalLength() { return cubic_spline_.GetTotalLength(); }

  py::dict getPolyline(const corridor::RealType delta_l) {
    // Constructs a Cartesian polyline from the cubic spline.
    // delta_s defines the sampling rate of the arc-length
    using namespace corridor;
    CartesianPoints2D polyline;
    RealType query_l = 0.0;
    RealType max_length = cubic_spline_.GetTotalLength();
    while (query_l <= max_length) {
      polyline.emplace_back(cubic_spline_.GetPositionAt(query_l));
      query_l += delta_l;
    }
    if (query_l > max_length) {
      // Add last point
      polyline.emplace_back(cubic_spline_.GetPositionAt(max_length));
    }
    // convert to python data structures
    const std::pair<py::list, py::list> py_polyline = to_py_lists(polyline);
    // Fill dictionary
    py::dict polyline_dict;
    polyline_dict["x"] = py_polyline.first;
    polyline_dict["y"] = py_polyline.second;
    return polyline_dict;
  }
};