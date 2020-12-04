#include "corridor/corridor_assignment/corridor_assignment.h"

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/corridor.h"
#include "corridor/frenet_types.h"
#include "corridor/internal/math.h"

namespace corridor {

CorridorRelatedFeatures ComputeCorridorRelatedObjectFeature(
    const CartesianState2D& cartesian_state,
    const BoxDimension& bounding_box_dimension, const Corridor& corridor) {
  CorridorRelatedFeatures features;
  features.frenet_frame = corridor.FrenetFrame(cartesian_state.position());
  features.frenet_state =
      features.frenet_frame.FromCartesianState(cartesian_state);

  // Direct access to the box dimensions (no transformation required)
  const auto projection_pair = bounding_box_dimension.projection(
      features.frenet_state.orientation().value);
  features.longitudinal_box_projection = projection_pair.first;
  features.lateral_box_projection = projection_pair.second;

  features.corridor_width =
      corridor.widthAt(features.frenet_frame.arc_length());

  features.corridor_length = corridor.lengthReferenceLine();
  features.corridor_center_offset =
      corridor.centerOffset(features.frenet_frame.arc_length());

  return features;
};

RealType ComputeAssignmentConfidence(const CorridorRelatedFeatures& features) {
  // Resulting likelihood that an object is located on the corridor as the
  // product of the likelihood of each feature.

  const auto latConf = LateralConfidence(features);
  const auto lonConf = 1.0;
  // LongitudinalConfidence(features);

  return latConf * lonConf;
};

// SemanticLabels ComputeMovingOnSemantics(
//     const CorridorRelatedFeatures& features) {
//   SemanticLabels longitudinal = LongitudinalMovingConfidence(features, 3.0);

//   SemanticLabels lateral = LateralMovingConfidence(features, 3.0);

//   return longitudinal + lateral;
// };

// //
// /////////////////////////////////////////////////////////////////////////////
// // Lateral and Longitudinal assignment confidences
// //
// /////////////////////////////////////////////////////////////////////////////

RealType LateralConfidence(const CorridorRelatedFeatures& features) {
  // Lateral features
  const RealType d = features.frenet_state.d();
  const RealType sigma_d =
      std::sqrt(features.frenet_state.covarianceMatrix().dd());

  const RealType object_width = features.lateral_box_projection;

  const RealType corridor_width = features.corridor_width;
  const RealType half_corridor_width = 0.5 * corridor_width;
  const RealType distance_from_center = d - features.corridor_center_offset;

  // std::cout << __FUNCTION__ << std::endl;
  // std::cout << "d = " << d << std::endl;
  // std::cout << "sigma_d = " << sigma_d << std::endl;
  // std::cout << "object_width = " << object_width << std::endl;
  // std::cout << "half_corridor_width = " << half_corridor_width <<
  // std::endl;
  // std::cout << "distance_from_center = " << distance_from_center <<
  // std::endl;

  RealType r1, r2, r3;
  // For the sake of optimization, we distinguish between the two cases:
  // Case 1): object's projection is almost like that of a point
  // Case 2): it's reasonably wide. In which case, we further distinguish
  // between two cases:
  // 		Case a): it's at least as wide as the corridor
  //		Case b): it's width is less than that of the corridor.
  if (object_width < std::numeric_limits<RealType>::epsilon()) {
    // Case 1) object projection is nearly 0: In this case, we should only
    // evaluate the middle integral, and avoid paying the cost of computing the
    // other two integrals:
    r1 = 0.0;
    r2 = math::evaluateIntegralLineWidthGaussian(0.0, 1.0, distance_from_center,
                                                 sigma_d, -half_corridor_width,
                                                 half_corridor_width);
    r3 = 0.0;
  } else {
    // Case 2) objects's projection is reasonably wide, in which case, we see if
    // the object's, perhaps, too wide:
    const RealType m = 1.0 / object_width;
    const RealType b = 0.5 * (1.0 + corridor_width / object_width);
    if (object_width >= corridor_width) {
      // Case 2.a): we only need to evaluate the first and the third integral
      r1 = math::evaluateIntegralLineWidthGaussian(
          m, b, distance_from_center, sigma_d,
          -0.5 * (corridor_width + object_width), 0.0);
      r2 = 0.0;
      r3 = math::evaluateIntegralLineWidthGaussian(
          -m, b, distance_from_center, sigma_d, 0.0,
          0.5 * (corridor_width + object_width));
    } else {
      // Case 2.b): we need to evaluate all three parts of the integral
      r1 = math::evaluateIntegralLineWidthGaussian(
          m, b, distance_from_center, sigma_d,
          -0.5 * (corridor_width + object_width),
          -0.5 * (corridor_width - object_width));
      r2 = math::evaluateIntegralLineWidthGaussian(
          0.0, 1.0, distance_from_center, sigma_d,
          -0.5 * (corridor_width - object_width),
          0.5 * (corridor_width - object_width));
      r3 = math::evaluateIntegralLineWidthGaussian(
          -m, b, distance_from_center, sigma_d,
          0.5 * (corridor_width - object_width),
          0.5 * (corridor_width + object_width));
    }
  }
  // std::cout << "r1, r2, r3 = " << r1 << ", " << r2 << ", " << r3 <<
  // std::endl;
  return r1 + r2 + r3;
};

RealType LongitudinalConfidence(const CorridorRelatedFeatures& features) {
  // Longitudinal features
  const RealType l = features.frenet_state.l();
  const RealType sigma_l =
      std::sqrt(features.frenet_state.covarianceMatrix().ll());
  const RealType object_length = features.longitudinal_box_projection;
  const RealType corridor_length = features.corridor_length;

  // std::cout << __FUNCTION__ << std::endl;
  // std::cout << "l = " << l << std::endl;
  // std::cout << "sigma_l = " << sigma_l << std::endl;
  // std::cout << "object_length = " << object_length << std::endl;
  // std::cout << "corridor_length = " << corridor_length << std::endl;

  // For the sake of optimization, we distinguish between the two cases:
  // Case 1): object's projection is almost like that of a point
  // Case 2): it's reasonably long.
  RealType r1, r2, r3;
  if (object_length <= std::numeric_limits<RealType>::epsilon()) {
    // Case 1) object projection is nearly 0: In this case, we should only
    // evaluate the middle integral, and avoid paying the cost of computing the
    // other two integrals:
    r1 = 0.0;
    r2 = math::evaluateIntegralLineWidthGaussian(0.0, 1.0, l, sigma_l, 0.0,
                                                 corridor_length);
    r3 = 0.0;
  } else {
    // Case 2) objects's projection is reasonably long, in which case, we see if
    // the object's, perhaps, too long:
    const RealType m1 = 1.0 / object_length;
    const RealType m2 = -m1;
    const RealType b1 = 0.5;
    const RealType b2 = 0.5 + corridor_length / object_length;

    const RealType half_obj_length = 0.5 * object_length;

    if (object_length >= corridor_length) {
      // Case 2.a): we only need to evaluate the first and the third integral
      r1 = math::evaluateIntegralLineWidthGaussian(
          m1, b1, l, sigma_l, -half_obj_length, 0.5 * corridor_length);
      r2 = 0.0;
      r3 = math::evaluateIntegralLineWidthGaussian(
          m2, b2, l, sigma_l, 0.5 * corridor_length,
          corridor_length + half_obj_length);
    } else {
      r1 = math::evaluateIntegralLineWidthGaussian(
          m1, b1, l, sigma_l, -half_obj_length, half_obj_length);
      r2 = math::evaluateIntegralLineWidthGaussian(
          0, 1, l, sigma_l, half_obj_length, corridor_length - half_obj_length);
      r3 = math::evaluateIntegralLineWidthGaussian(
          m2, b2, l, sigma_l, corridor_length - half_obj_length,
          corridor_length + half_obj_length);
    }
  }
  // std::cout << "r1, r2, r3 = " << r1 << ", " << r2 << ", " << r3 <<
  // std::endl;
  return r1 + r2 + r3;
};

// //
// /////////////////////////////////////////////////////////////////////////////
// // Lateral and Longitudinal moving confidences
// //
// /////////////////////////////////////////////////////////////////////////////

// SemanticLabels LongitudinalMovingConfidence(
//     const CorridorRelatedFeatures& features, const RealType sigma_band) {
//   // Longitudinal features
//   const RealType vl = features.frenet_state.vl();
//   const RealType sigma_vl =
//       std::sqrt(features.frenet_state_covMat.velocity().ll());

//   // Velocities below 0.5 mps are considered as non-moving
//   // TODO parameter!
//   const RealType nonMoving_limit = std::max(sigma_vl * sigma_band, 0.5);

//   SemanticLabels semantic_labels;
//   // Confidence based on negative velocity
//   semantic_labels[SemanticLabel::kUpstream] =
//       math::evaluateIntegralLineWidthGaussian(0.0, 1.0, vl, sigma_vl, -100.0,
//                                               -nonMoving_limit);

//   // Non moving confidence
//   semantic_labels[SemanticLabel::kDownstream] =
//       math::evaluateIntegralLineWidthGaussian(0.0, 1.0, vl, sigma_vl,
//                                               nonMoving_limit, 100.0);

//   // Confidence based on positive velocity
//   semantic_labels[SemanticLabel::kLongitudinalNonMoving] =
//       math::evaluateIntegralLineWidthGaussian(
//           0.0, 1.0, vl, sigma_vl, -nonMoving_limit, nonMoving_limit);

//   return semantic_labels;
// }

// SemanticLabels LateralMovingConfidence(const CorridorRelatedFeatures&
// features,
//                                        const RealType sigma_band) {
//   // Lateral features
//   const RealType vd = features.frenet_state.vd();
//   const RealType sigma_vd =
//       std::sqrt(features.frenet_state_covMat.velocity().dd());

//   // Velocities below 0.5 mps are considered as non-moving
//   // TODO parameter!
//   const RealType nonMoving_limit = std::max(sigma_vd * sigma_band, 0.5);

//   SemanticLabels semantic_labels;

//   // Confidence based on negative velocity
//   semantic_labels[SemanticLabel::kRight] =
//       math::evaluateIntegralLineWidthGaussian(0.0, 1.0, vd, sigma_vd, -100.0,
//                                               -nonMoving_limit);

//   // Non moving confidence based on negative velocity
//   semantic_labels[SemanticLabel::kLeft] =
//       math::evaluateIntegralLineWidthGaussian(0.0, 1.0, vd, sigma_vd,
//                                               nonMoving_limit, 100.0);

//   semantic_labels[SemanticLabel::kLateralNonMoving] =
//       math::evaluateIntegralLineWidthGaussian(
//           0.0, 1.0, vd, sigma_vd, -nonMoving_limit, nonMoving_limit);

//   return semantic_labels;
// }

}  // namespace corridor
