#pragma once

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/corridor.h"
#include "corridor/frenet_types.h"

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

  BoxDimension box_dimension;

  // TODO: use uncertainty value once the long and lat projection uncertainty is
  // needed (UT)
  RealType longitudinal_box_projection;
  RealType lateral_box_projection;

  // Corridor features
  RealType corridor_width;          //< at projection point
  RealType corridor_length;         //< from start to end
  RealType corridor_center_offset;  //< offset reference line to center

  bool longitudinalMatched() const {
    return (0.0 <= frenet_state.position().l()) &&
           (frenet_state.position().l() <= corridor_length);
  }

  CorridorRelatedFeatures(void)
      : frenet_frame(),
        frenet_state(),
        box_dimension(),
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
    const BoxDimension& bounding_box_dimension, const Corridor& corridor);

RealType LateralConfidence(const CorridorRelatedFeatures& features);

RealType LongitudinalConfidence(const CorridorRelatedFeatures& features);

RealType ComputeAssignmentConfidence(const CorridorRelatedFeatures& features);

}  // namespace corridor
