#include "corridor/frenet_types.h"

#include <Eigen/Core>
#include <iostream>

#include "corridor/cubic_spline/cubic_spline_types.h"
#include "corridor/unscented_transformation/polar_coordinate_transformation.h"

using namespace corridor;

// /////////////////////////////////////////////////////////////////////////////
// Frenet Frame
// /////////////////////////////////////////////////////////////////////////////

CartesianPoint2D FrenetFrame2D::FromFrenetPoint(
    const FrenetPoint2D& frenet_point) const {
  // Vector transfromation, translated by the origin of the frenet frame
  return FromFrenetVector(frenet_point) + origin_;
};

CartesianVector2D FrenetFrame2D::FromFrenetVector(
    const FrenetVector2D& frenet_vector) const {
  //! Local point
  const FrenetVector2D relative_vector{
      frenet_vector.l() - frenet_base_.arc_length, frenet_vector.d()};
  // Coordination transformation
  return rotMat_F2C_ * relative_vector;
};

FrenetPoint2D FrenetFrame2D::FromCartesianPoint(
    const CartesianPoint2D& cartesian_position) const {
  //! Define local point
  const CartesianVector2D relative_vector = cartesian_position - origin_;
  FrenetVector2D tmp = FromCartesianVector(relative_vector);
  // Add arc-length of the frenet base
  tmp.l() += frenet_base_.arc_length;
  //! Convert it to a FrenetPoint and return
  return FrenetPoint2D(tmp);
};

FrenetVector2D FrenetFrame2D::FromCartesianVector(
    const CartesianVector2D& cartesian_vector) const {
  //! Coordination transformation
  FrenetVector2D frenet_vector = rotMat_C2F_ * cartesian_vector;
  return frenet_vector;
};

RealType FrenetFrame2D::FromCartesianOrientation(
    const RealType cartesian_orientation) const {
  const auto angle = cartesian_orientation - frenet_base_.orientation;
  return constrainAngle(angle);
}

FrenetStateVector2D FrenetFrame2D::FromCartesianStateVector(
    const CartesianStateVector2D& cartesian_state,
    const bool moving_frenet_frame) const {
  // Simple transformation as projection of the cartesian state position and
  // vectors onto the axes of the Frenet frame
  FrenetPoint2D position = FromCartesianPoint(cartesian_state.position());
  FrenetVector2D velocity = FromCartesianVector(cartesian_state.velocity());

  if (moving_frenet_frame) {
    // In case of assuming a moving Frenet frame we need some correction terms
    // for accounting the rotation of the Frenet Frame.
    // Velocity and acceleration of the Frenet base is assumed to be the
    // projection of the velocity onto the tangent vector of the curve.
    const RealType vel_p = velocity.l();
    const RealType theta_dot = frenet_base_.curvature * vel_p;

    const CartesianVector2D relative_vector =
        cartesian_state.position() - origin_;
    RotationMatrix rotMat_C2F_prime;
    rotMat_C2F_prime << normal_.x(), normal_.y(), -tangent_.x(), -tangent_.y();
    velocity += theta_dot * rotMat_C2F_prime * relative_vector;
  }
  return FrenetStateVector2D(position, velocity);
}

FrenetStateCovarianceMatrix2D FrenetFrame2D::FromCartesianStateCovarianceMatrix(
    const CartesianStateCovarianceMatrix2D& state_vector_covariance_matrix,
    const bool moving_frenet_frame) const {
  // First simple quasi-linear transformation of the covariance matrices.
  const auto pos_covMat = state_vector_covariance_matrix.position();
  const auto vel_covMat = state_vector_covariance_matrix.velocity();
  const auto posvel_covMat = state_vector_covariance_matrix.pos_vel();

  const auto frenet_pos_covMat =
      rotMat_C2F_ * pos_covMat * rotMat_C2F_.transpose();
  const auto frenet_vel_covMat =
      rotMat_C2F_ * vel_covMat * rotMat_C2F_.transpose();
  const auto frenet_pos_vel_covMat =
      rotMat_C2F_ * posvel_covMat * rotMat_C2F_.transpose();

  return FrenetStateCovarianceMatrix2D(frenet_pos_covMat, frenet_vel_covMat,
                                       frenet_pos_vel_covMat);
};

FrenetState2D FrenetFrame2D::FromCartesianState(
    const CartesianState2D& cartesian_state, const bool moving_frenet_frame) {
  return {
      FromCartesianStateVector(cartesian_state.mean()),
      FromCartesianStateCovarianceMatrix(cartesian_state.covarianceMatrix())};
}

// /////////////////////////////////////////////////////////////////////////////
// Frenet Polyline
// /////////////////////////////////////////////////////////////////////////////

RealType FrenetPolyline::deviationAt(const RealType query_l) const {
  // Check if query l value is smaller or larger that the polyline
  const auto s_min = data_(DataType::kArclength, 0);
  const auto s_max = data_(DataType::kArclength, data_.cols() - 1);

  if (query_l <= (s_min + cubic_spline::g_epsilon_projection)) {
    return data_(DataType::kDeviation, 0);
  }
  if ((s_max - cubic_spline::g_epsilon_projection) <= query_l) {
    return data_(DataType::kDeviation, data_.cols() - 1);
  }

  //! Find index of closest data point to a given l value.
  auto start_index = 0;
  (data_.colwise() - FrenetPoint2D(query_l, 0))
      .topRows<1>()
      .colwise()
      .lpNorm<1>()
      .minCoeff(&start_index);

  //! If query_l is smaller then the l-value of the closest node, shift the
  //! start index by one
  const auto s_node = data_(DataType::kArclength, start_index);
  if (start_index > 0 && query_l < s_node) {
    start_index--;
  }

  if (start_index + 1 >= data_.cols()) {
    return data_(DataType::kDeviation, start_index);
  }

  const auto delta_l = data_(DataType::kArclength, start_index + 1) -
                       data_(DataType::kArclength, start_index);
  assert(abs(delta_l) > 1e-3);
  const auto alpha =
      (query_l - data_(DataType::kArclength, start_index)) / delta_l;

  const auto delta_d = data_(DataType::kDeviation, start_index + 1) -
                       data_(DataType::kDeviation, start_index);

  return delta_d * alpha + data_(DataType::kDeviation, start_index);
};

// /////////////////////////////////////////////////////////////////////////////
// Frenet state (mean and covariance matrix)
// /////////////////////////////////////////////////////////////////////////////

UncertainValue FrenetState2D::abs_velocity() {
  const auto& polar_velocity_state = polarVelocityState();
  return polar_velocity_state->abs_value();
};
UncertainValue FrenetState2D::orientation() {
  const auto& polar_velocity_state = polarVelocityState();
  return polar_velocity_state->orientation();
};

const PolarStatePtr FrenetState2D::polarVelocityState() {
  if (polar_velocity_state_ != nullptr) {
    // if polar velocity is already calculated return pointer.
    return polar_velocity_state_;
  }

  // Initialize shared ptr
  polar_velocity_state_ = std::make_shared<PolarState2D>();

  // Polar velocity is not yet set, calculate and return it.
  unscented_transformation::ToPolarCoordinates2D(
      mean_.velocity(), cov_mat_.velocity(), &polar_velocity_state_->mean,
      &polar_velocity_state_->cov_mat);

  return polar_velocity_state_;
};
