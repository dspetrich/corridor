#include "corridor/corridor_assignment/corridor_assignment.h"

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/corridor.h"
#include "corridor/corridor_assignment/corridor_related_semantics.h"
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

  const auto latConf = LateralAssignmentConfidence(features);
  const auto lonConf = LongitudinalAssignmentConfidence(features);

  return latConf * lonConf;
};

// /////////////////////////////////////////////////////////////////////////////
// Lateral and Longitudinal assignment confidences
// /////////////////////////////////////////////////////////////////////////////

RealType LateralAssignmentConfidence(const CorridorRelatedFeatures& features) {
  // Lateral features
  const RealType d = features.frenet_state.d();
  const RealType sigma_d =
      std::sqrt(features.frenet_state.covarianceMatrix().dd());

  const RealType object_width = features.lateral_box_projection;

  const RealType corridor_width = features.corridor_width;
  const RealType half_corridor_width = 0.5 * corridor_width;
  const RealType distance_from_center = d - features.corridor_center_offset;

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

RealType LongitudinalAssignmentConfidence(
    const CorridorRelatedFeatures& features) {
  // Longitudinal features
  const RealType l = features.frenet_state.l();
  const RealType sigma_l =
      std::sqrt(features.frenet_state.covarianceMatrix().ll());
  const RealType object_length = features.longitudinal_box_projection;
  const RealType corridor_length = features.corridor_length;

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

// /////////////////////////////////////////////////////////////////////////////
// Moving Confidence
// /////////////////////////////////////////////////////////////////////////////

RealType MovingConfidence(const UncertainValue& absolute_velocity,
                          const RealType sigma_band) {
  // Velocities below 0.1 mps are always considered as non-moving
  // TODO parameter!
  const RealType nonMoving_limit =
      std::max(absolute_velocity.standard_deviation * sigma_band, 0.1);

  // Moving confidence is 1-non_moving_confidence
  const RealType non_moving_confidence =
      math::evaluateIntegralLineWidthGaussian(
          0.0, 1.0, absolute_velocity.value,
          absolute_velocity.standard_deviation, -100, nonMoving_limit);

  // return 1.0 - non_moving_confidence;
  return 1 - non_moving_confidence;
};

// /////////////////////////////////////////////////////////////////////////////
// Relative Direction Confidence
// /////////////////////////////////////////////////////////////////////////////
SemanticLabelSet RelativeDirectionConfidence(
    const UncertainValue& relative_heading_angle, const RealType sigma_band) {
  // Extract value and standard deviation for better readability
  const RealType v = relative_heading_angle.value;
  const RealType sigma_v = relative_heading_angle.standard_deviation;

  // Deviation from straight forward/crossing angle which is still considered
  // forward/crossing
  // const static RealType delta_phi = M_PI / 64.0;  // ~ 2.8 degrees
  const static RealType delta_phi = std::max(sigma_v * sigma_band, M_PI / 64.0);

  // Specific angles
  const static RealType a0 = -10.0;
  const static RealType a1 = -2 * M_PI + delta_phi;
  const static RealType a2 = -1.5 * M_PI - delta_phi;
  const static RealType a3 = -1.5 * M_PI + delta_phi;
  const static RealType a4 = -M_PI - delta_phi;
  const static RealType a5 = -M_PI + delta_phi;
  const static RealType a6 = -0.5 * M_PI - delta_phi;
  const static RealType a7 = -0.5 * M_PI + delta_phi;
  const static RealType a8 = -delta_phi;
  const static RealType a9 = delta_phi;
  const static RealType a10 = -a7;
  const static RealType a11 = -a6;
  const static RealType a12 = -a5;
  const static RealType a13 = -a4;
  const static RealType a14 = -a3;
  const static RealType a15 = -a2;
  const static RealType a16 = -a1;
  const static RealType a17 = -a0;

  // Slope of the non-horizontal lines
  const static RealType denom = (M_PI - 4 * delta_phi);
  const static RealType m = 2.0 / denom;

  // b values for non-horizontal lines
  const static RealType b1 = 1.0;
  const static RealType b2 = -3 * (M_PI + 2. / 3. * delta_phi) / denom;
  const static RealType b3 = 2 * (2 * M_PI - delta_phi) / denom;
  const static RealType b4 = 1.0;
  const static RealType b5 = -2 * (M_PI + delta_phi) / denom;
  const static RealType b6 = 3 * (M_PI - 2. / 3. * delta_phi) / denom;
  const static RealType b7 = 1.0;
  const static RealType b8 = -(M_PI + 2 * delta_phi) / denom;
  const static RealType b9 = 2 * (M_PI - delta_phi) / denom;
  const static RealType b10 = 1.0;
  const static RealType b11 = -2 * delta_phi / denom;
  const static RealType b12 = (M_PI - 2 * delta_phi) / denom;
  const static RealType b13 = 1.0;
  const static RealType b14 = b12;
  const static RealType b15 = b11;
  const static RealType b16 = 1.0;
  const static RealType b17 = b9;
  const static RealType b18 = b8;
  const static RealType b19 = 1.0;
  const static RealType b20 = b6;
  const static RealType b21 = b5;
  const static RealType b22 = 1.0;
  const static RealType b23 = b3;
  const static RealType b24 = b2;
  const static RealType b25 = 1.0;

  // Initialize
  SemanticLabelSet semantic_labels(
      {SemanticLabel::kDownstream, SemanticLabel::kUpstream,
       SemanticLabel::kTowardsLeft, SemanticLabel::kTowardsRight});

  // A) following downstream
  RealType following_downstream = 0.0;
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(0.0, b1, v, sigma_v, a0, a1);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(-m, b2, v, sigma_v, a1, a2);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(m, b12, v, sigma_v, a7, a8);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(0.0, b13, v, sigma_v, a8, a9);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(-m, b14, v, sigma_v, a9, a10);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(m, b24, v, sigma_v, a15, a16);
  following_downstream +=
      math::evaluateIntegralLineWidthGaussian(0.0, b25, v, sigma_v, a16, a17);

  semantic_labels.setLabel(SemanticLabel::kDownstream, following_downstream);

  // B) following upstream
  RealType following_upstream = 0.0;
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(m, b6, v, sigma_v, a3, a4);
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(0.0, b7, v, sigma_v, a4, a5);
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(-m, b8, v, sigma_v, a5, a6);
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(m, b18, v, sigma_v, a11, a12);
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(0.0, b19, v, sigma_v, a12, a13);
  following_upstream +=
      math::evaluateIntegralLineWidthGaussian(-m, b20, v, sigma_v, a13, a14);

  semantic_labels.setLabel(SemanticLabel::kUpstream, following_upstream);

  // C) crossing towards left
  RealType crossing_left = 0.0;
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(m, b3, v, sigma_v, a1, a2);
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(0.0, b4, v, sigma_v, a2, a3);
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(-m, b5, v, sigma_v, a3, a4);
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(m, b15, v, sigma_v, a9, a10);
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(0.0, b16, v, sigma_v, a10, a11);
  crossing_left +=
      math::evaluateIntegralLineWidthGaussian(-m, b17, v, sigma_v, a11, a12);

  semantic_labels.setLabel(SemanticLabel::kTowardsLeft, crossing_left);

  // D) crossing towards right
  RealType crossing_right = 0.0;
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(m, b9, v, sigma_v, a5, a6);
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(0.0, b10, v, sigma_v, a6, a7);
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(-m, b11, v, sigma_v, a7, a8);
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(m, b21, v, sigma_v, a13, a14);
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(0.0, b22, v, sigma_v, a14, a15);
  crossing_right +=
      math::evaluateIntegralLineWidthGaussian(-m, b23, v, sigma_v, a15, a16);

  semantic_labels.setLabel(SemanticLabel::kTowardsRight, crossing_right);

  // Guarantee that the sum of all labels is one
  // semantic_labels.normalize();

  return semantic_labels;
}

}  // namespace corridor
