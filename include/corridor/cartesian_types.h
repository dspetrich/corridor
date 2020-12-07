#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iomanip>  // std::setprecision

#include "corridor/basic_types.h"

namespace corridor {

using CartesianPoint2D = BasicPoint2D;
using CartesianVector2D = BasicVector2D;

using CartesianPoints2D = std::vector<CartesianPoint2D>;
using CartesianVectors2D = std::vector<CartesianVector2D>;

// /////////////////////////////////////////////////////////////////////////////
// Cartesian state vector and covariances
// /////////////////////////////////////////////////////////////////////////////

struct CartesianStateVector2D
    : public Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign> {
  CartesianStateVector2D(void)
      : Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>() {}
  CartesianStateVector2D(const CartesianPoint2D& position,
                         const CartesianVector2D& velocity) {
    (*this) << position.x(), position.y(), velocity.x(), velocity.y();
  }
  CartesianStateVector2D(const RealType x, const RealType y,
                         const RealType vx = 0.0, const RealType vy = 0.0) {
    (*this) << x, y, vx, vy;
  }

  using Base = Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>;

  // This constructor allows you to construct CartesianStateVector2D from Eigen
  // expressions
  template <typename OtherDerived>
  CartesianStateVector2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to
  // CartesianStateVector2D
  template <typename OtherDerived>
  CartesianStateVector2D& operator=(
      const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Copy
  CartesianPoint2D position() const { return this->head<2>(); }
  CartesianVector2D velocity() const { return this->tail<2>(); }

  // Non-mutable views
  const RealType x() const { return (*this)[0]; }
  const RealType y() const { return (*this)[1]; }
  const RealType vx() const { return (*this)[2]; }
  const RealType vy() const { return (*this)[3]; }

  // Mutable views
  RealType& x() { return (*this)[0]; }
  RealType& y() { return (*this)[1]; }
  RealType& vx() { return (*this)[2]; }
  RealType& vy() { return (*this)[3]; }
};

struct CovarianceMatrix2D
    : public Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign> {
  CovarianceMatrix2D(const RealType xx = 1e-12, const RealType yy = 1e-12,
                     const RealType xy = 0.0) {
    (*this) << xx, xy, xy, yy;
  }

  using Base = Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>;

  // This constructor allows you to construct CovarianceMatrix2D from Eigen
  // expressions
  template <typename OtherDerived>
  CovarianceMatrix2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to CovarianceMatrix2D
  template <typename OtherDerived>
  CovarianceMatrix2D& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Non-mutable views
  const RealType xx() const { return (*this)(0, 0); }
  const RealType yy() const { return (*this)(1, 1); }
  const RealType xy() const { return (*this)(0, 1); }
};

struct CartesianStateCovarianceMatrix2D
    : public Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign> {
  CartesianStateCovarianceMatrix2D(void)
      : CartesianStateCovarianceMatrix2D(CovarianceMatrix2D(),
                                         CovarianceMatrix2D()) {}
  CartesianStateCovarianceMatrix2D(
      const CovarianceMatrix2D& pos, const CovarianceMatrix2D& vel,
      const CovarianceMatrix2D& pos_vel = CovarianceMatrix2D::Zero()) {
    this->block<2, 2>(0, 0) = pos;
    this->block<2, 2>(2, 2) = vel;
    this->block<2, 2>(0, 2) = pos_vel;
    this->block<2, 2>(2, 0) = pos_vel;
  }

  using Base = Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign>;

  // This constructor allows you to construct CartesianStateCovarianceMatrix2D
  // from Eigen expressions
  template <typename OtherDerived>
  CartesianStateCovarianceMatrix2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to
  // CartesianStateCovarianceMatrix2D
  template <typename OtherDerived>
  CartesianStateCovarianceMatrix2D& operator=(
      const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Copies
  CovarianceMatrix2D position() const { return this->block<2, 2>(0, 0); }
  CovarianceMatrix2D velocity() const { return this->block<2, 2>(2, 2); }
  CovarianceMatrix2D pos_vel() const { return this->block<2, 2>(0, 2); }
};

// /////////////////////////////////////////////////////////////////////////////
// Cartesian state (mean and covariance matrix)
// /////////////////////////////////////////////////////////////////////////////

class CartesianState2D {
 public:
  CartesianState2D(void) : mean_(), cov_mat_() {}
  CartesianState2D(
      const CartesianPoint2D& position, const CartesianVector2D& velocity,
      const CovarianceMatrix2D& cm_position = CovarianceMatrix2D(),
      const CovarianceMatrix2D& cm_velocity = CovarianceMatrix2D(),
      const CovarianceMatrix2D& cm_pos_vel = CovarianceMatrix2D::Zero())
      : mean_(position, velocity),
        cov_mat_(cm_position, cm_velocity, cm_pos_vel) {}

  // Simple getter
  RealType x() const { return mean_.x(); }
  RealType y() const { return mean_.y(); }
  RealType vx() const { return mean_.vx(); }
  RealType vy() const { return mean_.vy(); }

  CartesianPoint2D position() const { return mean_.position(); }
  CartesianPoint2D velocity() const { return mean_.velocity(); }
  const CartesianStateVector2D& mean() const { return mean_; }
  const CartesianStateCovarianceMatrix2D& covarianceMatrix() const {
    return cov_mat_;
  }

  UncertainValue abs_velocity() {
    const auto& polar_velocity_state = getPolarVelocityStatePtr();
    return polar_velocity_state->abs_value();
  }
  UncertainValue orientation() {
    const auto& polar_velocity_state = getPolarVelocityStatePtr();
    return polar_velocity_state->orientation();
  }

 private:
  CartesianStateVector2D mean_;
  CartesianStateCovarianceMatrix2D cov_mat_;

  // optional polar interpretation of the velocity vector. Will only be
  // constructed if and then cached.
  const PolarStatePtr getPolarVelocityStatePtr();
  PolarStatePtr polar_velocity_state_;
};

/**
 * @brief Simple 2D object class.
 * Mostly for example purposes of the corridor assignment function
 *
 */
class CartesianObjectState2D {
 public:
  CartesianObjectState2D(void)
      : id_(InvalidId), box_dimension_(), state_(), state_cov_mat_() {}

  CartesianPoint2D centerPosition() const { return state_.position(); }
  const BoxDimension& dimension() const { return box_dimension_; }

  const CartesianStateVector2D& state() const { return state_; }
  const CartesianStateCovarianceMatrix2D& stateCovarianceMatrix() const {
    return state_cov_mat_;
  }

  UncertainValue occupancyRadius() const {
    using namespace std;
    // Radius of the box, with the uncertainty of the position
    UncertainValue radius = box_dimension_.enclosingCircleRadius();
    const RealType position_variance = state_cov_mat_.position().trace();
    radius.standard_deviation += std::sqrt(position_variance);
    return radius;
  }

  /**
   * @brief Returns min and max point of the bounding box, multiplied by the
   * sigma band of the radius uncertainty (includes position uncertainty)
   *
   * @param sigma_band: sigma of the radius dimension uncertainty
   * @return std::pair<CartesianPoint2D, CartesianPoint2D>
   */

  AlignedBoundingBox2D boundingBox(const RealType sigma_band = 1e-12) const {
    const UncertainValue radius = occupancyRadius();
    const RealType delta_value =
        M_SQRT1_2 * (radius.value + sigma_band * radius.standard_deviation);
    const CartesianPoint2D delta(delta_value, delta_value);
    AlignedBoundingBox2D abb;
    abb.min() = state_.position() - delta;
    abb.max() = state_.position() + delta;
    return abb;
  }

 private:
  IdType id_;
  BoxDimension box_dimension_;
  CartesianStateVector2D state_;
  CartesianStateCovarianceMatrix2D state_cov_mat_;
};

}  // namespace corridor