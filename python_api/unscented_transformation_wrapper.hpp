#pragma once

#include <boost/python.hpp>

#include "corridor/unscented_transformation/polar_coordinate_transformation.h"
#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

// /////////////////////////////////////////////////////////////////////////////
// Uncertainty Transformation
// /////////////////////////////////////////////////////////////////////////////

namespace py = boost::python;

struct FlatCartesianStateAndCovMat2D {
  corridor::RealType x;
  corridor::RealType y;
  corridor::RealType var_x;
  corridor::RealType var_y;
  corridor::RealType cov_xy;
};

struct FlatPolarStateAndCovMat2D {
  corridor::RealType r;
  corridor::RealType phi;
  corridor::RealType var_r;
  corridor::RealType var_phi;
  corridor::RealType cov_rphi;
};

corridor::StateMeanAndCovarianceMatrix Convert(
    const FlatCartesianStateAndCovMat2D& flat_cartesian_state) {
  corridor::StateMeanAndCovarianceMatrix state(2);
  state.mean << flat_cartesian_state.x, flat_cartesian_state.y;
  state.covMat(0, 0) = flat_cartesian_state.var_x;
  state.covMat(1, 1) = flat_cartesian_state.var_y;
  state.covMat(1, 0) = flat_cartesian_state.cov_xy;
  state.covMat(0, 1) = flat_cartesian_state.cov_xy;
  return state;
}

FlatPolarStateAndCovMat2D Convert(
    const corridor::StateMeanAndCovarianceMatrix& polar_state) {
  FlatPolarStateAndCovMat2D flat_polar_state;
  flat_polar_state.r = polar_state.mean(0);
  flat_polar_state.phi = polar_state.mean(1);
  flat_polar_state.var_r = polar_state.covMat(0, 0);
  flat_polar_state.var_phi = polar_state.covMat(1, 1);
  flat_polar_state.cov_rphi = polar_state.covMat(1, 0);
  return flat_polar_state;
};

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

FlatPolarStateAndCovMat2D UnscentedTransformationPolarCoordinate2D(
    const FlatCartesianStateAndCovMat2D& cartesian_state) {
  using namespace corridor;
  StateMeanAndCovarianceMatrix polar_state =
      UnscentedTransformationPolarCoordinates(Convert(cartesian_state));
  return Convert(polar_state);
}
