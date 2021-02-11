#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <boost/optional.hpp>
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

// To be replaced by std::optional once c++17 is established
template <typename T>
using Optional = boost::optional<T>;

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

// TODO: template!
struct UncertainValue {
  UncertainValue(const RealType _value = 0.0,
                 const RealType _standard_deviation = 1e-12)
      : value(_value), standard_deviation(_standard_deviation) {}

  RealType variance() const { return standard_deviation * standard_deviation; }
  RealType value;  // mean value
  RealType standard_deviation;
};
// introspection
inline std::ostream &operator<<(std::ostream &os, const UncertainValue &uv) {
  os << "value = " << uv.value << "; variance = " << uv.variance();
  return os;
};

// TODO: maybe add more functionality later and make it an orientated bounding
// box if needed.
using AlignedBoundingBox2D = Eigen::AlignedBox<RealType, 2>;

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
      : BasicPoint2D(abs_value, constrainAngle(orientation)) {}

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
  // marginalization of absolute value
  UncertainValue abs_value() const {
    return {mean.abs_value(), std::sqrt(cov_mat.var_abs_value())};
  }
  // marginalization of orientation angle
  UncertainValue orientation() const {
    return {mean.orientation(), std::sqrt(cov_mat.var_orientation())};
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