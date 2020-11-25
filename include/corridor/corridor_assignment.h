#pragma once

// #include "situation_analysis/SceneObject.h"
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
  FrenetStateVector2D frenet_state;
  FrenetStateCovarianceMatrix2D frenet_state_covMat;
  SceneObjectBoxShape frenet_shape;
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
        frenet_state_covMat(),
        frenet_shape(),
        corridor_width(-1),
        corridor_length(-1),
        corridor_center_offset(-1) {}
};
// introspection
inline std::ostream& operator<<(std::ostream& os,
                                const CorridorRelatedFeatures& crf) {
  os << "Corridor Assignment Features\n";
  os << crf.frenet_frame;
  os << "Frenet State: " << crf.frenet_state.transpose() << "\n";
  os << crf.frenet_shape;
  os << "corridor_width: " << crf.corridor_width << "\n";
  os << "corridor_length: " << crf.corridor_length << "\n";
  os << "corridor_center_offset: " << crf.corridor_center_offset << "\n";
  return os;
};

CorridorRelatedFeatures ComputeCorridorRelatedObjectFeature(
    const SceneObject& object, const Corridor& corridor);

RealType LateralConfidence(const CorridorRelatedFeatures& features);

RealType LongitudinalConfidence(const CorridorRelatedFeatures& features);

RealType ComputeAssignmentConfidence(const CorridorRelatedFeatures& features);

}  // namespace corridor
