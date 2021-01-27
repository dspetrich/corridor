#pragma once

#include <boost/python.hpp>

#include "corridor/corridor.h"

// Phython API
#include "cubic_spline_wrapper.hpp"
#include "utility.hpp"

namespace py = boost::python;

// /////////////////////////////////////////////////////////////////////////////
// Corridor wrapper
// /////////////////////////////////////////////////////////////////////////////

struct CorridorWrapper {
  // CPP object
  corridor::Corridor corridor_;

  CorridorWrapper(const int id, const py::list& node_x,
                  const py::list& node_y) {
    using namespace corridor;
    CartesianPoints2D points = toCartesianPoints(node_x, node_y);
    const RealType distance_left_boundary = 2.0;
    const RealType distance_right_boundary = 2.0;
    corridor_ =
        Corridor(id, points, distance_left_boundary, distance_right_boundary);
  }

  CorridorWrapper(const int id, const py::list& node_x, const py::list& node_y,
                  const py::list& py_first_tangent,
                  const py::list& py_last_tangent) {
    using namespace corridor;
    const CartesianPoints2D points = toCartesianPoints(node_x, node_y);
    const CartesianVector2D first_tangent =
        to_cartesian_vector(py_first_tangent);
    const CartesianVector2D last_tangent = to_cartesian_vector(py_last_tangent);
    const RealType distance_left_boundary = 2.0;
    const RealType distance_right_boundary = 2.0;
    corridor_ = Corridor(id, points, first_tangent, last_tangent,
                         distance_left_boundary, distance_right_boundary);
  }

  py::dict GetCartesianPolylinesLines(const corridor::RealType delta_s) {
    using namespace corridor;
    // Construct Cartesian polylines form corridors reference line and
    // boundaries. delta_s defines the sampling rate for all three lines
    CartesianPoints2D reference_line;
    CartesianPoints2D left_boundary;
    CartesianPoints2D right_boundary;

    corridor_.fillCartesianPolylines(delta_s, &reference_line, &left_boundary,
                                     &right_boundary);

    // convert to python data structures
    const std::pair<py::list, py::list> py_refline =
        to_py_lists(reference_line);
    const std::pair<py::list, py::list> py_left_boundary =
        to_py_lists(left_boundary);
    const std::pair<py::list, py::list> py_right_boundary =
        to_py_lists(right_boundary);

    // Fill dictionary
    py::dict polylines;
    polylines["reference_line_x"] = py_refline.first;
    polylines["reference_line_y"] = py_refline.second;
    polylines["left_boundary_x"] = py_left_boundary.first;
    polylines["left_boundary_y"] = py_left_boundary.second;
    polylines["right_boundary_x"] = py_right_boundary.first;
    polylines["right_boundary_y"] = py_right_boundary.second;
    return polylines;
  }

  py::list ToFrenetStateVector(const py::list& py_cartesian_state_vector,
                               const bool moving_frenet_frame) {
    using namespace corridor;
    CartesianStateVector2D cartesian_state_vector =
        convert(py_cartesian_state_vector);

    FrenetFrame2D frenet_frame =
        corridor_.FrenetFrame(cartesian_state_vector.position());

    FrenetStateVector2D frenet_state = frenet_frame.FromCartesianStateVector(
        cartesian_state_vector, moving_frenet_frame);

    return convert(frenet_state);
  }

  FlatFrenetStateAndCovMat2D ToFrenetState(
      const FlatCartesianStateAndCovMat2D& flat_cartesian_state,
      const bool moving_frenet_frame) {
    using namespace corridor;
    CartesianState2D cartesian_state = Convert(flat_cartesian_state);

    FrenetFrame2D frenet_frame =
        corridor_.FrenetFrame(cartesian_state.position());

    FrenetState2D frenet_state =
        frenet_frame.FromCartesianState(cartesian_state, moving_frenet_frame);

    return Convert(frenet_state);
  }

  py::dict GetFrenetFrame(const cr::RealType target_x,
                          const cr::RealType target_y) {
    using namespace corridor;
    CartesianPoint2D target_point(target_x, target_y);

    FrenetFrame2D frenet_frame = corridor_.FrenetFrame(target_point);

    return to_py_dict(frenet_frame);
  }

  corridor::RealType lengthReferenceLine() {
    return corridor_.lengthReferenceLine();
  }

  corridor::RealType curvatureAt(const corridor::RealType l) {
    return corridor_.curvatureAt(l);
  }
};

void TestCorridorHandle(const CorridorWrapper& corridor_wrapper) {
  std::cout << corridor_wrapper.corridor_.id() << std::endl;
}
