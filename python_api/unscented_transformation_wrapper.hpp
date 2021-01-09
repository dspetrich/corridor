#pragma once

#include <boost/python.hpp>

#include "corridor/unscented_transformation/polar_coordinate_transformation.h"
#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/state_transformation.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

// Phython API
#include "corridor_wrapper.hpp"
#include "utility.hpp"

namespace py = boost::python;

// /////////////////////////////////////////////////////////////////////////////
// Polar Coordinate transformation
// /////////////////////////////////////////////////////////////////////////////

corridor::StateMeanAndCovarianceMatrix UnscentedTransformationPolarCoordinates(
    const corridor::StateMeanAndCovarianceMatrix& initial_state) {
  namespace ut = corridor::unscented_transformation;
  corridor::PolarVector2D polar_mean;
  corridor::PolarCovarianceMatrix2D polar_cov_mat;

  ut::ToPolarCoordinates2D(initial_state.mean, initial_state.covMat,
                           &polar_mean, &polar_cov_mat);

  return corridor::StateMeanAndCovarianceMatrix(polar_mean, polar_cov_mat);
};

py::list pyCartesianToPolarTransformation2D(const corridor::RealType x1,
                                            const corridor::RealType x2) {
  namespace ut = corridor::unscented_transformation;
  Eigen::Vector2d polar_vector = ut::CartesianToPolarTransformation2D({x1, x2});
  py::list polar_coordinates;
  polar_coordinates.append(polar_vector(0));
  polar_coordinates.append(polar_vector(1));
  return polar_coordinates;
}

py::list pyPolarToCartesianTransformation2D(const corridor::RealType radius,
                                            const corridor::RealType phi) {
  namespace ut = corridor::unscented_transformation;
  Eigen::Vector2d cartesian_vector =
      ut::PolarToCartesianTransformation2D({radius, phi});
  py::list cartesian_coordinates;
  cartesian_coordinates.append(cartesian_vector(0));
  cartesian_coordinates.append(cartesian_vector(1));
  return cartesian_coordinates;
}

FlatPolarPositionAndCovMat2D UnscentedTransformationPolarCoordinate2D(
    const FlatCartesianPositionAndCovMat2D& cartesian_state) {
  using namespace corridor;
  StateMeanAndCovarianceMatrix polar_state =
      UnscentedTransformationPolarCoordinates(Convert(cartesian_state));
  return Convert(polar_state);
}

// /////////////////////////////////////////////////////////////////////////////
// State transformation
// /////////////////////////////////////////////////////////////////////////////

FlatFrenetStateAndCovMat2D UnscentedStateTransformation(
    const CorridorWrapper& corridor_wrapper,
    const FlatCartesianStateAndCovMat2D flat_cartesian_state,
    const bool moving_frenet_frame) {
  using namespace corridor;
  namespace ut = unscented_transformation;
  FrenetState2D frenet_state = ut::ToFrenetState(corridor_wrapper.corridor_,
                                                 Convert(flat_cartesian_state));
  return Convert(frenet_state);
}