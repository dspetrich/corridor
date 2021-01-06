#pragma once

#include <boost/python.hpp>

#include "corridor/corridor.h"
#include "cubic_spline_wrapper.hpp"

// /////////////////////////////////////////////////////////////////////////////
// Utility functions
// /////////////////////////////////////////////////////////////////////////////

// Create Corridor from list of points

struct CorridorWrapper {
  corridor::Corridor corridor;

  CorridorWrapper(const int id, const py::list& node_x,
                  const py::list& node_y) {
    using namespace corridor;
    CartesianPoints2D points = toCartesianPoints(node_x, node_y);
    const RealType distance_left_boundary = 2.0;
    const RealType distance_right_boundary = 2.0;
    corridor =
        Corridor(id, points, distance_left_boundary, distance_right_boundary);
  }

  void ExtraxtCartesianLines(const corridor::RealType delta_s) {
    // Extract spline coefficients from corridor, convert left and right
    // boundary at delta_s
  }
};

void TestCorridorHandle(const CorridorWrapper& corridor_wrapper) {
  std::cout << corridor_wrapper.corridor.id() << std::endl;
}
