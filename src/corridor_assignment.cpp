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
    const OrientedBoundingBox& oriented_bounding_box,
    const Corridor& corridor) {
  CorridorRelatedFeatures features;
  features.frenet_frame = corridor.FrenetFrame(cartesian_state.position());
  features.frenet_state =
      features.frenet_frame.FromCartesianState(cartesian_state);

  // Direct access to the box dimensions (no transformation required)
  const UncertainValue orientation = oriented_bounding_box.orientation();
  const auto relative_orientation =
      orientation.value - features.frenet_frame.frenet_base().orientation;
  features.frenet_obb = OrientedBoundingBox(
      {relative_orientation, orientation.standard_deviation},
      oriented_bounding_box.length(), oriented_bounding_box.width());

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

  const RealType object_width = features.frenet_obb.projectionX2();

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
    // evaluate the middle integral, and avoid paying the cost of computing
    // the other two integrals:
    r1 = 0.0;
    r2 = math::evaluateIntegralLineWidthGaussian(0.0, 1.0, distance_from_center,
                                                 sigma_d, -half_corridor_width,
                                                 half_corridor_width);
    r3 = 0.0;
  } else {
    // Case 2) objects's projection is reasonably wide, in which case, we see
    // if the object's, perhaps, too wide:
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

// Assignment confidence that the object is left of the corridor
RealType LeftLateralAssignmentConfidence(
    const CorridorRelatedFeatures& features) {
  // Lateral features
  const RealType d = features.frenet_state.d();
  const RealType sigma_d =
      std::sqrt(features.frenet_state.covarianceMatrix().dd());

  const RealType object_width = features.frenet_obb.projectionX2();

  const RealType corridor_width = features.corridor_width;
  const RealType half_corridor_width = 0.5 * corridor_width;
  const RealType distance_from_center = d - features.corridor_center_offset;

  RealType r1, r2;
  // For the sake of optimization, we distinguish between the two cases:
  // Case 1): object's projection is almost like that of a point
  // Case 2): it's reasonably wide. In which case, we further distinguish
  // between two cases:
  // 		Case a): it's at least as wide as the corridor
  //		Case b): it's width is less than that of the corridor.
  if (object_width < std::numeric_limits<RealType>::epsilon()) {
    // Case 1) object projection is nearly 0: In this case, we should only
    // evaluate the middle integral, and avoid paying the cost of computing
    // the other two integrals:

    r1 = math::evaluateIntegralLineWidthGaussian(
        0.0, 1.0, distance_from_center, sigma_d, -100, -half_corridor_width);
    r2 = 0.0;
  } else {
    // Case 2) objects's projection is reasonably wide, in which case, we see
    // if the object's, perhaps, too wide:
    const RealType m = -1.0 / object_width;
    const RealType b = (object_width - corridor_width) / (2.0 * object_width);

    // Case 2.b): we need to evaluate all three parts of the integral
    r1 = math::evaluateIntegralLineWidthGaussian(
        0.0, 1.0, distance_from_center, sigma_d, -100,
        -0.5 * (corridor_width + object_width));
    r2 = math::evaluateIntegralLineWidthGaussian(
        m, b, distance_from_center, sigma_d,
        -0.5 * (corridor_width + object_width),
        -0.5 * (corridor_width - object_width));
  }
  // std::cout << "r1, r2, r3 = " << r1 << ", " << r2 << ", " << r3 <<
  // std::endl;
  return r1 + r2;
}

// Assignment confidence that the object is right of the corridor
RealType RightLateralAssignmentConfidence(
    const CorridorRelatedFeatures& features) {
  // Lateral features
  const RealType d = features.frenet_state.d();
  const RealType sigma_d =
      std::sqrt(features.frenet_state.covarianceMatrix().dd());

  const RealType object_width = features.frenet_obb.projectionX2();

  const RealType corridor_width = features.corridor_width;
  const RealType half_corridor_width = 0.5 * corridor_width;
  const RealType distance_from_center = d - features.corridor_center_offset;

  RealType r1, r2;
  // For the sake of optimization, we distinguish between the two cases:
  // Case 1): object's projection is almost like that of a point
  // Case 2): it's reasonably wide. In which case, we further distinguish
  // between two cases:
  // 		Case a): it's at least as wide as the corridor
  //		Case b): it's width is less than that of the corridor.
  if (object_width < std::numeric_limits<RealType>::epsilon()) {
    // Case 1) object projection is nearly 0: In this case, we should only
    // evaluate the middle integral, and avoid paying the cost of computing
    // the other two integrals:
    r1 = 0.0;
    r2 = math::evaluateIntegralLineWidthGaussian(
        0.0, 1.0, distance_from_center, sigma_d, half_corridor_width, 100.0);

  } else {
    // Case 2) objects's projection is reasonably wide, in which case, we see
    // if the object's, perhaps, too wide:
    const RealType m = 1.0 / object_width;
    const RealType b = (object_width - corridor_width) / (2.0 * object_width);

    // Case 2.b): we need to evaluate all three parts of the integral
    r1 = math::evaluateIntegralLineWidthGaussian(
        m, b, distance_from_center, sigma_d,
        0.5 * (corridor_width - object_width),
        0.5 * (corridor_width + object_width));
    r2 = math::evaluateIntegralLineWidthGaussian(
        0.0, 1.0, distance_from_center, sigma_d,
        0.5 * (corridor_width + object_width), 100.0);
  }
  // std::cout << "r1, r2, r3 = " << r1 << ", " << r2 << ", " << r3 <<
  // std::endl;
  return r1 + r2;
}

RealType LongitudinalAssignmentConfidence(
    const CorridorRelatedFeatures& features) {
  // Longitudinal features
  const RealType l = features.frenet_state.l();
  const RealType sigma_l =
      std::sqrt(features.frenet_state.covarianceMatrix().ll());
  const RealType object_length = features.frenet_obb.projectionX1();
  const RealType corridor_length = features.corridor_length;

  // For the sake of optimization, we distinguish between the two cases:
  // Case 1): object's projection is almost like that of a point
  // Case 2): it's reasonably long.
  RealType r1, r2, r3;
  if (object_length <= std::numeric_limits<RealType>::epsilon()) {
    // Case 1) object projection is nearly 0: In this case, we should only
    // evaluate the middle integral, and avoid paying the cost of computing
    // the other two integrals:
    r1 = 0.0;
    r2 = math::evaluateIntegralLineWidthGaussian(0.0, 1.0, l, sigma_l, 0.0,
                                                 corridor_length);
    r3 = 0.0;
  } else {
    // Case 2) objects's projection is reasonably long, in which case, we see
    // if the object's, perhaps, too long:
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

  // Moving confidence
  const RealType moving_confidence = math::evaluateIntegralLineWidthGaussian(
      0.0, 1.0, absolute_velocity.value, absolute_velocity.standard_deviation,
      nonMoving_limit, 100);

  return moving_confidence;
};

// /////////////////////////////////////////////////////////////////////////////
// Relative Direction Confidence
// /////////////////////////////////////////////////////////////////////////////

RealType RelativeOrientationConfidence(
    const RealType direction_angle,
    const UncertainValue& relative_heading_angle, const RealType sigma_band) {
  // Extract value and standard deviation for better readability
  const RealType relative_phi =
      constrainAngle(relative_heading_angle.value - direction_angle);
  const RealType sigma_phi = relative_heading_angle.standard_deviation;

  // Deviation from straight forward/crossing angle which is still considered
  // forward/crossing
  // const static RealType delta_phi = M_PI / 64.0;  // ~ 2.8 degrees
  const static RealType delta_phi =
      std::max(sigma_phi * sigma_band, M_PI / 64.0);

  const RealType denominator = (M_PI - 4 * delta_phi);
  const RealType m = 2.0 / denominator;
  const RealType b = (M_PI - 2 * delta_phi) / denominator;

  const RealType p1 = -M_PI_2 + delta_phi;
  const RealType p2 = -delta_phi;
  const RealType p3 = delta_phi;
  const RealType p4 = M_PI_2 - delta_phi;

  const RealType r1 = math::evaluateIntegralLineWidthGaussian(
      m, b, relative_phi, sigma_phi, p1, p2);
  const RealType r2 = math::evaluateIntegralLineWidthGaussian(
      0, 1, relative_phi, sigma_phi, p2, p3);
  const RealType r3 = math::evaluateIntegralLineWidthGaussian(
      -m, b, relative_phi, sigma_phi, p3, p4);

  return r1 + r2 + r3;
}

}  // namespace corridor
