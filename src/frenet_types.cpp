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

FrenetState2D FrenetFrame2D::FromCartesianState(
    const CartesianState2D& cartesian_state,
    const bool moving_frenet_frame) const {
  if (moving_frenet_frame) {
    return FromCartesianStateTaylorExpansion(cartesian_state);
  }
  // Linear transformation function
  return {
      this->frenet_base().id, FromCartesianStateVector(cartesian_state.mean()),
      FromCartesianStateCovarianceMatrix(cartesian_state.covarianceMatrix())};
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
    // Velocity of Frenet base is assumed to be the projection of the velocity
    // onto the tangent vector of the curve.
    const RealType theta_dot = frenet_base_.curvature * velocity.l();
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
  // 4x4 rotation matrix
  Eigen::Matrix<corridor::RealType, 4, 4> rotation_matrix =
      Eigen::Matrix<corridor::RealType, 4, 4>::Zero();

  // Fill with 2d rotation matrices
  rotation_matrix.block<2, 2>(0, 0) = rotMat_C2F_;
  rotation_matrix.block<2, 2>(2, 2) = rotMat_C2F_;

  // Linear transformation of the covariance matrices.
  return rotation_matrix * state_vector_covariance_matrix *
         rotation_matrix.transpose();
}

FrenetState2D FrenetFrame2D::FromCartesianStateTaylorExpansion(
    const CartesianState2D& cartesian_state) const {
  // Non-linear transformation using the Taylor Series upt to term
  FrenetStateVector2D frenet_state_vector =
      FromCartesianStateVector(cartesian_state.mean(), true);

  JacobianMatrix jacobian_matrix = defineJacobianMatrix(cartesian_state);
  FrenetStateCovarianceMatrix2D frenet_cov_mat =
      jacobian_matrix * cartesian_state.covarianceMatrix() *
      jacobian_matrix.transpose();

  return {this->frenet_base().id, frenet_state_vector, frenet_cov_mat};
}

FrenetFrame2D::JacobianMatrix FrenetFrame2D::defineJacobianMatrix(
    const CartesianState2D& cartesian_state) const {
  // Easy access
  const CartesianVector2D relative_vector =
      cartesian_state.position() - origin_;
  const RealType projection_on_tangent = tangent_.dot(relative_vector);
  const RealType projection_on_normal = normal_.dot(relative_vector);

  const RealType vp = tangent_.dot(cartesian_state.velocity());

  // Jacobean matrix at cartesian mean
  Eigen::Matrix<RealType, 4, 4> jacobian_matrix =
      Eigen::Matrix<RealType, 4, 4>::Zero();
  jacobian_matrix.block<2, 2>(0, 0) = rotMat_C2F_;
  jacobian_matrix(2, 0) = frenet_base_.curvature * vp * normal_.x();
  jacobian_matrix(2, 1) = frenet_base_.curvature * vp * normal_.y();
  jacobian_matrix(2, 2) = tangent_.x() + frenet_base_.curvature * tangent_.x() *
                                             projection_on_normal;
  jacobian_matrix(2, 3) = tangent_.y() + frenet_base_.curvature * tangent_.y() *
                                             projection_on_normal;

  jacobian_matrix(3, 0) = -frenet_base_.curvature * vp * tangent_.x();
  jacobian_matrix(3, 1) = -frenet_base_.curvature * vp * tangent_.y();
  jacobian_matrix(3, 2) = normal_.x() + frenet_base_.curvature * tangent_.x() *
                                            projection_on_tangent;
  jacobian_matrix(3, 3) = normal_.y() + frenet_base_.curvature * tangent_.y() *
                                            projection_on_tangent;

  return jacobian_matrix;
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

  //! Get index of segment which contains the arc-length
  DataMatrix::Index index = 0;
  DataMatrix::Index max_index =
      data_.cols() - 2;  // first index of last segment
  const bool valid = (data_.row(kArclength).array() > query_l).maxCoeff(&index);

  if (valid) {
    // Check that always a valid data segment is used
    index = (0 < index) ? (index - 1) : (index);
    index = (max_index < index) ? (max_index) : (index);
  } else {
    // Arc-length is longer as the refernce line. Use last segment start
    // index
    index = (valid) ? (index) : (max_index);
  }

  if (index + 1 >= data_.cols()) {
    return data_(DataType::kDeviation, index);
  }

  // Interpolation of the deviation
  const auto delta_l = data_(DataType::kArclength, index + 1) -
                       data_(DataType::kArclength, index);
  assert(abs(delta_l) > 1e-3);
  const auto alpha = (query_l - data_(DataType::kArclength, index)) / delta_l;

  const auto delta_d = data_(DataType::kDeviation, index + 1) -
                       data_(DataType::kDeviation, index);

  return delta_d * alpha + data_(DataType::kDeviation, index);
};

// /////////////////////////////////////////////////////////////////////////////
// Frenet state (mean and covariance matrix)
// /////////////////////////////////////////////////////////////////////////////

UncertainValue FrenetState2D::abs_velocity() {
  const auto& polar_velocity_state = getPolarVelocityStatePtr();
  return polar_velocity_state->abs_value();
};
UncertainValue FrenetState2D::orientation() {
  const auto& polar_velocity_state = getPolarVelocityStatePtr();
  return polar_velocity_state->orientation();
};

const PolarStatePtr FrenetState2D::getPolarVelocityStatePtr() {
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
