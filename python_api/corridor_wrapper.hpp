#pragma once

#include <boost/python.hpp>

#include "corridor/corridor.h"
#include "cubic_spline_wrapper.hpp"

namespace py = boost::python;

// /////////////////////////////////////////////////////////////////////////////
// Utility functions
// /////////////////////////////////////////////////////////////////////////////

corridor::CartesianStateVector2D convert(
    const py::list& cartesian_state_vector) {
  using namespace corridor;
  return CartesianStateVector2D(
      py::extract<RealType>(cartesian_state_vector[0]),
      py::extract<RealType>(cartesian_state_vector[1]),
      py::extract<RealType>(cartesian_state_vector[2]),
      py::extract<RealType>(cartesian_state_vector[3]));
}

py::list convert(const corridor::FrenetStateVector2D& frenet_state_vector) {
  py::list py_frenet_state_vector;
  py_frenet_state_vector.append(frenet_state_vector.l());
  py_frenet_state_vector.append(frenet_state_vector.d());
  py_frenet_state_vector.append(frenet_state_vector.vl());
  py_frenet_state_vector.append(frenet_state_vector.vd());
  return py_frenet_state_vector;
}

struct CorridorWrapper {
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

  py::list ToFrenetStateVector(const py::list& py_cartesian_state_vector) {
    using namespace corridor;
    CartesianStateVector2D cartesian_state_vector =
        convert(py_cartesian_state_vector);

    FrenetFrame2D frenet_frame =
        corridor_.FrenetFrame(cartesian_state_vector.position());

    FrenetStateVector2D frenet_state =
        frenet_frame.FromCartesianStateVector(cartesian_state_vector);

    return convert(frenet_state);
  }
};

void TestCorridorHandle(const CorridorWrapper& corridor_wrapper) {
  std::cout << corridor_wrapper.corridor_.id() << std::endl;
}
