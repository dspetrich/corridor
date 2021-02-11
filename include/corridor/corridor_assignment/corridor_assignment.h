#pragma once

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/corridor.h"
#include "corridor/corridor_assignment/corridor_related_semantics.h"
#include "corridor/frenet_types.h"
#include "corridor/internal/oriented_bounding_box.h"

namespace corridor {

/**
 * @brief Collection of object features computed relatively to one corridor
 * structure. The information is given in the Frenet frame of the reference line
 * at the projection's point of the object's reference point.
 *
 */
struct CorridorRelatedFeatures {
  // Object features
  FrenetFrame2D frenet_frame;
  FrenetState2D frenet_state;

  // TODO: use uncertainty value once the long and lat projection uncertainty is
  // needed (UT)
  RealType longitudinal_box_projection;
  RealType lateral_box_projection;

  // Corridor features
  RealType corridor_width;          //< at projection point
  RealType corridor_length;         //< from start to end
  RealType corridor_center_offset;  //< offset reference line to center

  bool longitudinallyMatched() const {
    // TODO include uncertainty
    return (0.0 <= frenet_state.position().l()) &&
           (frenet_state.position().l() <= corridor_length);
  }

  UncertainValue longitudinalVelocity() const {
    return {frenet_state.vl(), frenet_state.covarianceMatrix().vlvl()};
  }

  UncertainValue lateralVelocity() const {
    return {frenet_state.vd(), frenet_state.covarianceMatrix().vdvd()};
  }

  CorridorRelatedFeatures(void)
      : frenet_frame(),
        frenet_state(),
        corridor_width(-1),
        corridor_length(-1),
        corridor_center_offset(-1) {}
};
// introspection
inline std::ostream& operator<<(std::ostream& os,
                                const CorridorRelatedFeatures& crf) {
  os << "Corridor Assignment Features\n";
  os << crf.frenet_frame;
  os << "Frenet State: " << crf.frenet_state.mean().transpose() << "\n";
  os << "box projections (lat,long): " << crf.lateral_box_projection << ", "
     << crf.longitudinal_box_projection << "\n";
  os << "corridor_width: " << crf.corridor_width << "\n";
  os << "corridor_length: " << crf.corridor_length << "\n";
  os << "corridor_center_offset: " << crf.corridor_center_offset << "\n";
  return os;
};

CorridorRelatedFeatures ComputeCorridorRelatedObjectFeature(
    const CartesianState2D& cartesian_state,
    const OrientedBoundingBox& oriented_bounding_box, const Corridor& corridor);

// /////////////////////////////////////////////////////////////////////////////
// // Lateral and Longitudinal assignment confidences
// /////////////////////////////////////////////////////////////////////////////

RealType LateralAssignmentConfidence(const CorridorRelatedFeatures& features);

RealType LongitudinalAssignmentConfidence(
    const CorridorRelatedFeatures& features);

RealType ComputeAssignmentConfidence(const CorridorRelatedFeatures& features);

// /////////////////////////////////////////////////////////////////////////////
// Moving Confidences
// /////////////////////////////////////////////////////////////////////////////

/**
 * @brief Overall moving confidence of the object, based on the objects absolute
 * velocity and velocity uncertainty. If the absolute velocity is small compared
 * to the uncertainty, the moving confidence is small as well.
 *
 * @param cartesian_state: current cartesian motion state with covariance matrix
 * @return RealType: moving confidence between 0 and 1.
 */
RealType MovingConfidence(const UncertainValue& absolute_velocity,
                          const RealType sigma_band);

RealType RelativeOrientationConfidence(
    const RealType direction_angle,
    const UncertainValue& relative_heading_angle, const RealType sigma_band);

SemanticLabelSet RelativeDirectionConfidence(
    const UncertainValue& relative_heading_angle, const RealType sigma_band);

}  // namespace corridor
