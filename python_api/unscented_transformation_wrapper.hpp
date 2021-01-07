#pragma once

#include <boost/python.hpp>

#include "corridor/unscented_transformation/polar_coordinate_transformation.h"
#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/state_transformation.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

// Phython API
#include "corridor_wrapper.hpp"

namespace py = boost::python;

// /////////////////////////////////////////////////////////////////////////////
// Utility functions
// /////////////////////////////////////////////////////////////////////////////

struct FlatCartesianPositionAndCovMat2D {
  corridor::RealType x;
  corridor::RealType y;
  corridor::RealType var_x;
  corridor::RealType var_y;
  corridor::RealType cov_xy;
};

struct FlatCartesianStateAndCovMat2D {
  corridor::RealType x;
  corridor::RealType y;
  corridor::RealType vx;
  corridor::RealType vy;

  corridor::RealType var_x;
  corridor::RealType var_y;
  corridor::RealType var_vx;
  corridor::RealType var_vy;

  corridor::RealType cov_xy;
  corridor::RealType cov_xvx;
  corridor::RealType cov_xvy;
  corridor::RealType cov_yvx;
  corridor::RealType cov_yvy;
  corridor::RealType cov_vxvy;
};

struct FlatPolarPositionAndCovMat2D {
  corridor::RealType r;
  corridor::RealType phi;
  corridor::RealType var_r;
  corridor::RealType var_phi;
  corridor::RealType cov_rphi;
};

struct FlatFrenetStateAndCovMat2D {
  corridor::RealType l;
  corridor::RealType d;
  corridor::RealType vl;
  corridor::RealType vd;

  corridor::RealType var_l;
  corridor::RealType var_d;
  corridor::RealType var_vl;
  corridor::RealType var_vd;

  corridor::RealType cov_ld;
  corridor::RealType cov_lvl;
  corridor::RealType cov_lvd;
  corridor::RealType cov_dvl;
  corridor::RealType cov_dvd;
  corridor::RealType cov_vlvd;
};

corridor::StateMeanAndCovarianceMatrix Convert(
    const FlatCartesianPositionAndCovMat2D& flat_cartesian_state) {
  corridor::StateMeanAndCovarianceMatrix state(2);
  state.mean << flat_cartesian_state.x, flat_cartesian_state.y;
  state.covMat(0, 0) = flat_cartesian_state.var_x;
  state.covMat(1, 1) = flat_cartesian_state.var_y;
  state.covMat(1, 0) = flat_cartesian_state.cov_xy;
  state.covMat(0, 1) = flat_cartesian_state.cov_xy;
  return state;
}

corridor::CartesianState2D Convert(
    const FlatCartesianStateAndCovMat2D& flat_cartesian_state) {
  using namespace corridor;
  const CartesianPoint2D position(flat_cartesian_state.x,
                                  flat_cartesian_state.y);
  const CartesianVector2D velocity(flat_cartesian_state.vx,
                                   flat_cartesian_state.vy);
  const CovarianceMatrix2D cm_position(flat_cartesian_state.var_x,
                                       flat_cartesian_state.var_y,
                                       flat_cartesian_state.cov_xy);
  const CovarianceMatrix2D cm_velocity(flat_cartesian_state.var_vx,
                                       flat_cartesian_state.var_vy,
                                       flat_cartesian_state.cov_vxvy);
  const CovarianceMatrix2D cm_pos_vel(
      flat_cartesian_state.cov_xvx, flat_cartesian_state.cov_yvy,
      flat_cartesian_state.cov_xvy, flat_cartesian_state.cov_yvx);
  return CartesianState2D(position, velocity, cm_position, cm_velocity,
                          cm_pos_vel);
}

FlatFrenetStateAndCovMat2D Convert(
    const corridor::FrenetState2D& frenet_state) {
  FlatFrenetStateAndCovMat2D flat_frenet_state;
  flat_frenet_state.l = frenet_state.l();
  flat_frenet_state.d = frenet_state.d();
  flat_frenet_state.vl = frenet_state.vl();
  flat_frenet_state.vd = frenet_state.vd();

  const corridor::FrenetStateCovarianceMatrix2D& cov_mat =
      frenet_state.covarianceMatrix();

  flat_frenet_state.var_l = cov_mat.ll();
  flat_frenet_state.var_d = cov_mat.dd();
  flat_frenet_state.cov_ld = cov_mat.ld();

  flat_frenet_state.var_vl = cov_mat.vlvl();
  flat_frenet_state.var_vd = cov_mat.vdvd();
  flat_frenet_state.cov_vlvd = cov_mat.vlvd();

  flat_frenet_state.cov_lvl = cov_mat.lvl();
  flat_frenet_state.cov_lvd = cov_mat.lvd();
  flat_frenet_state.cov_dvl = cov_mat.dvl();
  flat_frenet_state.cov_dvd = cov_mat.dvd();

  return flat_frenet_state;
}

FlatPolarPositionAndCovMat2D Convert(
    const corridor::StateMeanAndCovarianceMatrix& polar_state) {
  FlatPolarPositionAndCovMat2D flat_polar_state;
  flat_polar_state.r = polar_state.mean(0);
  flat_polar_state.phi = polar_state.mean(1);
  flat_polar_state.var_r = polar_state.covMat(0, 0);
  flat_polar_state.var_phi = polar_state.covMat(1, 1);
  flat_polar_state.cov_rphi = polar_state.covMat(1, 0);
  return flat_polar_state;
};

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
    const FlatCartesianStateAndCovMat2D flat_cartesian_state) {
  using namespace corridor;
  namespace ut = unscented_transformation;
  FrenetState2D frenet_state = ut::ToFrenetState(corridor_wrapper.corridor_,
                                                 Convert(flat_cartesian_state));
  return Convert(frenet_state);
}