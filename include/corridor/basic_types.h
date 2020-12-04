#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace corridor {

// /////////////////////////////////////////////////////////////////////////////
// Basic types
// /////////////////////////////////////////////////////////////////////////////

using IdType = int64_t;
constexpr IdType InvalidId = -1;
using Ids = std::vector<IdType>;

using IdxType = size_t;

using RealType = double;

// Currently double, lets find a better representation later
using TimeStampType = double;
constexpr TimeStampType InvalidTimeStamp = -1;

/**
 * @brief Limits the value x to the interval [u,o] or [o,u], depending on which
 *        of both limits is greater.
 *
 * @tparam T: typename
 * @param x: value
 * @param u: first limit
 * @param o: second limit
 * @return T: constrained value
 */
template <typename T>
inline T limit(const T x, const T u, const T o) {
  return (o >= u) ? ((x <= o) ? ((x >= u) ? x : u) : o)
                  : ((x <= u) ? ((x >= o) ? x : o) : u);
};

/**
 * @brief Constraints the provided angle between -pi and pi
 *
 * @param angle: unconstrained angle
 * @return RealType: constrained angle
 */
inline RealType constrainAngle(RealType angle) {
  static const RealType two_pi = 2 * M_PI;
  angle = std::fmod(angle + M_PI, two_pi);
  if (angle < 0.0) {
    angle += two_pi;
  }
  return angle - M_PI;
}

// Can be easily templated
struct UncertainValue {
  UncertainValue(const RealType _value = 0.0, const RealType _variance = 1e-12)
      : value(_value), variance(_variance) {}

  RealType standardDeviation() const { return std::sqrt(variance); }
  RealType value;     // mean value
  RealType variance;  // variance = standardDeviation^2
};
// introspection
inline std::ostream &operator<<(std::ostream &os, const UncertainValue &uv) {
  os << "value = " << uv.value << "; variance = " << uv.variance << "\n";
  return os;
};

// TODO: maybe add more functionality later and make it an orientated bounding
// box if needed.
using AlignedBoundingBox2D = Eigen::AlignedBox<RealType, 2>;

/**
 * @brief Shape information of an object, as a bounding box as approximation
 * of its shape.
 *
 */
struct BoxDimension {
  UncertainValue length;
  UncertainValue width;

  BoxDimension(const RealType _length = 0.0, const RealType _width = 0.0)
      : length(_length), width(_width) {}

  BoxDimension(const UncertainValue &_length, const UncertainValue &_width)
      : length(_length), width(_width) {}

  BoxDimension(const BoxDimension &other)
      : length(other.length), width(other.width) {}

  /**
   * @brief Projection of the shape onto the x1 and x2 axis in which the
   * orientation angle is given. It is assumed that the orientation is defined
   * mathematically with respect to the x1 axis.
   *
   * @return std::pair<RealType, RealType>: first: projection onto x1 axis
   *                                        second: projection onto x2 axis
   */
  std::pair<RealType, RealType> projection(const RealType orientation) const {
    const auto cos = std::cos(orientation);
    const auto sin = std::sin(orientation);
    Eigen::Matrix2d projectionMatrix;
    projectionMatrix << cos, sin, sin, cos;
    const Eigen::Vector2d dimension(length.value, width.value);
    const Eigen::Vector2d projection = projectionMatrix * dimension;
    return std::make_pair(std::abs(projection.x()), std::abs(projection.y()));
  }

  RealType projectionX1(const RealType orientation) const {
    const Eigen::Vector2d projection(std::cos(orientation),
                                     std::sin(orientation));
    const Eigen::Vector2d dimension(length.value, width.value);
    return std::abs(projection.dot(dimension));
  }

  RealType projectionX2(const RealType orientation) const {
    const Eigen::Vector2d projection(std::sin(orientation),
                                     std::cos(orientation));
    const Eigen::Vector2d dimension(length.value, width.value);
    return std::abs(projection.dot(dimension));
  }

  UncertainValue enclosingCircleRadius() const noexcept {
    const RealType squared_length = length.value * length.value;
    const RealType squared_width = width.value * width.value;
    const RealType squared_value = 0.25 * (squared_length + squared_width);

    // Assuming length and width are uncorrelated
    const RealType variance =
        0.25 * ((squared_length / squared_value) * length.variance +
                (squared_width / squared_value) * width.variance);

    return {std::sqrt(squared_value), variance};
  }
};
// introspection
inline std::ostream &operator<<(std::ostream &os, const BoxDimension &box) {
  os << "\tlength: " << box.length;
  os << "\twidth: " << box.width << "\n";
  return os;
};

// /////////////////////////////////////////////////////////////////////////////
// Basic 2D point and vector
// /////////////////////////////////////////////////////////////////////////////

struct BasicPoint2D : public Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> {
  BasicPoint2D(void) : Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign>() {
    this->fill(0.0);
  }
  BasicPoint2D(const RealType x, const RealType y) { (*this) << x, y; }

  typedef Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> Base;

  // This constructor allows you to construct BasicPoint2D from Eigen
  // expressions
  template <typename OtherDerived>
  BasicPoint2D(const Eigen::MatrixBase<OtherDerived> &other)
      : Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to CartesianPoint2D
  template <typename OtherDerived>
  BasicPoint2D &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Base::operator=(other);
    return *this;
  }
};

/**
 * @brief Vector is a direction, point a position (with respect to a origin)
 */
using BasicVector2D = BasicPoint2D;

// /////////////////////////////////////////////////////////////////////////////
// Uncertainty state representation with mean and covariance matrix
// /////////////////////////////////////////////////////////////////////////////

// TODO: template this structure and the whole UT functions
struct StateMeanAndCovarianceMatrix {
  StateMeanAndCovarianceMatrix(const int state_dim = 0)
      : mean(Eigen::VectorXd(state_dim)),
        covMat(Eigen::MatrixXd(state_dim, state_dim)) {
    mean.fill(0.0);
    covMat.fill(0.0);
  }
  StateMeanAndCovarianceMatrix(const Eigen::VectorXd &_mean,
                               const Eigen::MatrixXd &_covMat)
      : mean(_mean), covMat(_covMat) {}
  Eigen::VectorXd mean;
  Eigen::MatrixXd covMat;
};

// Introspection
inline std::ostream &operator<<(std::ostream &os,
                                const StateMeanAndCovarianceMatrix &state) {
  os << "StateMeanAndCovarianceMatrix\n";
  os << "mean: " << state.mean.transpose() << "\n";
  os << state.covMat;
  return os;
}

// /////////////////////////////////////////////////////////////////////////////
// Polar 2D coordinates and covariance matrix
// /////////////////////////////////////////////////////////////////////////////

struct PolarVector2D : public BasicPoint2D {
  PolarVector2D(void) : BasicPoint2D() {}
  PolarVector2D(const RealType abs_value, const RealType orientation)
      : BasicPoint2D(abs_value, orientation) {}

  // This constructor allows you to construct FrenetPoint2D from Eigen
  // expressions
  template <typename OtherDerived>
  PolarVector2D(const Eigen::MatrixBase<OtherDerived> &other)
      : BasicPoint2D(other) {}

  typedef Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> Base;

  // This method allows you to assign Eigen expressions to FrenetPoint2D
  template <typename OtherDerived>
  PolarVector2D &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Base::operator=(other);
    return *this;
  }

  // Easy access to the polar coordinates
  const RealType abs_value() const { return (*this)[0]; }
  const RealType orientation() const { return (*this)[1]; }

  RealType &abs_value() { return (*this)[0]; }
  RealType &orientation() { return (*this)[1]; }
};

struct PolarCovarianceMatrix2D
    : public Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign> {
  PolarCovarianceMatrix2D(const RealType var_abs_value = 1e-12,
                          const RealType var_orientation = 1e-12,
                          const RealType cov_value_orientation = 0.0) {
    (*this) << var_abs_value, cov_value_orientation, cov_value_orientation,
        var_orientation;
  }

  using Base = Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>;

  // This constructor allows you to construct CovarianceMatrix2D from Eigen
  // expressions
  template <typename OtherDerived>
  PolarCovarianceMatrix2D(const Eigen::MatrixBase<OtherDerived> &other)
      : Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to CovarianceMatrix2D
  template <typename OtherDerived>
  PolarCovarianceMatrix2D &operator=(
      const Eigen::MatrixBase<OtherDerived> &other) {
    this->Base::operator=(other);
    return *this;
  }

  // Non-mutable views
  const RealType var_abs_value() const { return (*this)(0, 0); }
  const RealType var_orientation() const { return (*this)(1, 1); }
  const RealType cov_value_orientation() const { return (*this)(0, 1); }
};

struct PolarState2D {
  PolarVector2D mean;
  PolarCovarianceMatrix2D cov_mat;

  // Easy access functions:
  // a) marginalization
  UncertainValue abs_value() const {
    return {mean.abs_value(), cov_mat.var_abs_value()};
  }
  UncertainValue orientation() const {
    return {mean.orientation(), cov_mat.var_orientation()};
  }
};

using PolarStatePtr = std::shared_ptr<PolarState2D>;

// Introspection
inline std::ostream &operator<<(std::ostream &os, const PolarStatePtr state) {
  using namespace std;
  os << "Polar State Vector (|r|, phi): ";
  os << state->mean.transpose() << "\n";
  os << "Polar State CovMat:\n";
  os << state->cov_mat << "\n";
  return os;
};

}  // namespace corridor
