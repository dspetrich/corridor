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

struct CubicSpline {
  py::dict naturalSplineParameter(const py::list& node_x,
                                  const py::list& node_y) {
    using namespace corridor;

    CartesianPoints2D points = toCartesianPoints(node_x, node_y);

    // Create natural spline from points
    data_matrix = cs::NaturalSplineDataMatrixFromPoints(points);

    // Extract arc length information
    py::list arc_lengths;
    for (cs::DataIdx i = 0, max_idx = data_matrix.cols(); i < max_idx; i++) {
      arc_lengths.append(data_matrix.col(i)[cs::BasicDataTypes::kArcLength]);
    }

    // Construct spline coefficients
    cs::SplineCoefficients2d coefficients =
        cs::SplineCoefficientsFromDataMatrix(data_matrix);

    py::dict py_dict = to_py_dict(coefficients);
    py_dict["arc_length"] = arc_lengths;

    return py_dict;
  }

  py::dict clampedSplineParameter(const py::list& node_x,
                                  const py::list& node_y,
                                  const py::list& py_first_tangent,
                                  const py::list& py_last_tangent) {
    using namespace corridor;
    CartesianPoints2D points = toCartesianPoints(node_x, node_y);

    CartesianVector2D first_tangent = to_cartesian_vector(py_first_tangent);
    CartesianVector2D last_tangent = to_cartesian_vector(py_last_tangent);

    // Create clamped spline from points and tangents
    data_matrix = cs::ClampedSplineDataMatrixFromPoints(points, first_tangent,
                                                        last_tangent);

    // Extract arc length information
    py::list arc_lengths;
    for (cs::DataIdx i = 0, max_idx = data_matrix.cols(); i < max_idx; i++) {
      arc_lengths.append(data_matrix.col(i)[cs::BasicDataTypes::kArcLength]);
    }

    // Construct spline coefficients
    cs::SplineCoefficients2d coefficients =
        cs::SplineCoefficientsFromDataMatrix(data_matrix);

    py::dict py_dict = to_py_dict(coefficients);
    py_dict["arc_length"] = arc_lengths;

    return py_dict;
  }

  py::list constructFrenetFrames(const cr::RealType target_x,
                                 const cr::RealType target_y) {
    using namespace corridor;

    if (data_matrix.cols() == 0) {
      std::cout << "Data matrix not yet set! Call 'CreateSplineParams' first."
                << std::endl;
      return py::list();
    }

    CartesianPoint2D target_point;
    target_point << target_x, target_y;

    const FrenetFrames2D frenet_frames =
        cs::ConstructFrenetFrames(data_matrix, target_point);

    return to_py_list(frenet_frames);
  }

  // Data storage for augmented support points
  cs::DataMatrix<cr::RealType> data_matrix;
};