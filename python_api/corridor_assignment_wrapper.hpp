#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "corridor/basic_types.h"
#include "corridor/corridor.h"
#include "corridor/corridor_assignment/corridor_assignment.h"
#include "corridor/corridor_assignment/corridor_related_semantics.h"
#include "corridor/internal/math.h"

namespace py = boost::python;

// /////////////////////////////////////////////////////////////////////////////
// Utilities
// /////////////////////////////////////////////////////////////////////////////

struct FlatCorridorRelatedFeatures {
  corridor::RealType l;         //< arc-length of the object (mutable)
  corridor::RealType d;         //< deviation from the reference line (mutable)
  corridor::RealType sigma_l;   //< standard deviation of the arc-length
  corridor::RealType sigma_d;   //< standard deviation of the deviation
  corridor::RealType vl;        //< arc-length of the object (mutable)
  corridor::RealType vd;        //< deviation from the reference line (mutable)
  corridor::RealType sigma_vl;  //< standard deviation of the arc-length
  corridor::RealType sigma_vd;  //< standard deviation of the deviation
  corridor::RealType obj_width_ratio;   //< from the absolute corridor width
  corridor::RealType obj_length_ratio;  //< from the absolute corridor length

  // Corridor features
  corridor::RealType corridor_width;   //< at projection point
  corridor::RealType corridor_length;  //< from start to end
};

corridor::CorridorRelatedFeatures From(
    const FlatCorridorRelatedFeatures& flat_features) {
  using namespace corridor;
  CorridorRelatedFeatures features;

  const FrenetStateVector2D state_vector(flat_features.l, flat_features.d,
                                         flat_features.vl, flat_features.vd);
  const FrenetStateCovarianceMatrix2D state_cov_mat(
      flat_features.sigma_l * flat_features.sigma_l,
      flat_features.sigma_d * flat_features.sigma_d,
      flat_features.sigma_vl * flat_features.sigma_vl,
      flat_features.sigma_vd * flat_features.sigma_vd);

  features.frenet_state = FrenetState2D(state_vector, state_cov_mat);

  features.corridor_width = flat_features.corridor_width;
  features.corridor_length = flat_features.corridor_length;
  features.corridor_center_offset = 0.0;

  features.longitudinal_box_projection =
      (flat_features.corridor_length * flat_features.obj_length_ratio);
  features.lateral_box_projection =
      (flat_features.corridor_width * flat_features.obj_width_ratio);

  return features;
};

corridor::RealType LateralAssignmentConfidence(
    const FlatCorridorRelatedFeatures& flat_features) {
  return corridor::LateralAssignmentConfidence(From(flat_features));
};

corridor::RealType LongitudinalAssignmentConfidence(
    const FlatCorridorRelatedFeatures& flat_features) {
  return corridor::LongitudinalAssignmentConfidence(From(flat_features));
};

corridor::RealType evaluateIntegralLineWidthGaussian(
    const corridor::RealType m, const corridor::RealType b,
    const corridor::RealType x, const corridor::RealType sigma_original,
    const corridor::RealType lower_bound,
    const corridor::RealType upper_bound) {
  return corridor::math::evaluateIntegralLineWidthGaussian(
      m, b, x, sigma_original, lower_bound, upper_bound);
};

corridor::RealType MovingConfidence(const corridor::RealType abs_velocity,
                                    const corridor::RealType std_velocity,
                                    const corridor::RealType sigma_band) {
  return corridor::MovingConfidence({abs_velocity, std_velocity}, sigma_band);
}

py::dict RelativeDirectionConfidence(const corridor::RealType heading_mean,
                                     const corridor::RealType heading_std,
                                     const corridor::RealType sigma_band) {
  using namespace corridor;

  py::dict result;
  result["downstream"] = RelativeOrientationConfidence(
      0.0, {heading_mean, heading_std}, sigma_band);
  result["upstream"] = RelativeOrientationConfidence(
      M_PI, {heading_mean, heading_std}, sigma_band);
  result["towards_left"] = RelativeOrientationConfidence(
      M_PI_2, {heading_mean, heading_std}, sigma_band);
  result["towards_right"] = RelativeOrientationConfidence(
      -M_PI_2, {heading_mean, heading_std}, sigma_band);

  return result;
}