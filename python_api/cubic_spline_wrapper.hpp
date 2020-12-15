#pragma once

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numpy.hpp>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace cr = corridor;
namespace cs = cr::cubic_spline;

// /////////////////////////////////////////////////////////////////////////////
// Utility functions
// /////////////////////////////////////////////////////////////////////////////

template <typename T>
inline std::vector<T> to_std_vector(const py::object& iterable) {
  return std::vector<T>(py::stl_input_iterator<T>(iterable),
                        py::stl_input_iterator<T>());
}

inline cr::CartesianVector2D to_cartesian_vector(const py::list& py_list) {
  using namespace corridor;
  std::vector<RealType> std_vec = to_std_vector<RealType>(py_list);
  return CartesianVector2D(std_vec.front(), std_vec.back());
}

inline py::dict to_py_dict(const cs::SplineCoefficients2d& coeffs) {
  py::dict dict;

  py::list a_x, b_x, c_x, d_x;
  py::list a_y, b_y, c_y, d_y;
  for (const auto& c : coeffs) {
    a_x.append(c.a_x);
    b_x.append(c.b_x);
    c_x.append(c.c_x);
    d_x.append(c.d_x);

    a_y.append(c.a_y);
    b_y.append(c.b_y);
    c_y.append(c.c_y);
    d_y.append(c.d_y);
  }
  dict["a_x"] = a_x;
  dict["b_x"] = b_x;
  dict["c_x"] = c_x;
  dict["d_x"] = d_x;

  dict["a_y"] = a_y;
  dict["b_y"] = b_y;
  dict["c_y"] = c_y;
  dict["d_y"] = d_y;

  return dict;
}

inline py::dict to_py_dict(const cr::CartesianPoint2D& vector2d) {
  py::dict py_dict;
  py_dict["x"] = vector2d.x();
  py_dict["y"] = vector2d.y();
  return py_dict;
}

inline py::list to_py_list(const cr::FrenetFrames2D& frenet_frames) {
  py::list py_frenet_frames;
  for (const auto& ff : frenet_frames) {
    py::dict py_frenet_frame;
    py_frenet_frame["segm_index"] = ff.frenet_base().segment_info.idx;
    py_frenet_frame["segm_arc_length"] =
        ff.frenet_base().segment_info.relative_arc_length;
    py_frenet_frame["origin"] = to_py_dict(ff.origin());
    py_frenet_frame["tangent"] = to_py_dict(ff.tangent());
    py_frenet_frame["normal"] = to_py_dict(ff.normal());
    py_frenet_frames.append(py_frenet_frame);
  }
  return py_frenet_frames;
}

// /////////////////////////////////////////////////////////////////////////////
// Cubic Spline wrapper
// /////////////////////////////////////////////////////////////////////////////

struct CubicSpline {
  py::dict naturalSplineParameter(const py::list& node_x,
                                  const py::list& node_y) {
    using namespace corridor;
    std::vector<RealType> x_vec = to_std_vector<RealType>(node_x);
    std::vector<RealType> y_vec = to_std_vector<RealType>(node_y);

    cr::CartesianPoints2D points;
    for (int i = 0, n = x_vec.size(); i < n; i++) {
      points.emplace_back(x_vec[i], y_vec[i]);
    }
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

  py::dict clampedSplineParameter(const py::list& py_node_x,
                                  const py::list& py_node_y,
                                  const py::list& py_first_tangent,
                                  const py::list& py_last_tangent) {
    using namespace corridor;
    std::vector<RealType> x_vec = to_std_vector<RealType>(py_node_x);
    std::vector<RealType> y_vec = to_std_vector<RealType>(py_node_y);
    CartesianPoints2D points;
    for (int i = 0, n = x_vec.size(); i < n; i++) {
      points.emplace_back(x_vec[i], y_vec[i]);
    }

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