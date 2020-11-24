#pragma once

#include <Eigen/Core>
#include <iomanip>  // std::setprecision

#include "corridor/basic_types.h"

namespace corridor {

struct BasicPoint2D : public Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> {
  BasicPoint2D(void) : Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign>() {
    this->fill(0.0);
  }
  BasicPoint2D(const RealType x, const RealType y) { (*this) << x, y; }

  typedef Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> Base;

  // This constructor allows you to construct BasicPoint2D from Eigen
  // expressions
  template <typename OtherDerived>
  BasicPoint2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to CartesianPoint2D
  template <typename OtherDerived>
  BasicPoint2D& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
};

/**
 * @brief Vector is a direction, point a position (with respect to a origin)
 */
using BasicVector2D = BasicPoint2D;

using CartesianPoint2D = BasicPoint2D;
using CartesianVector2D = BasicVector2D;

using CartesianPoints2D = std::vector<CartesianPoint2D>;
using CartesianVectors2D = std::vector<CartesianVector2D>;

}  // namespace corridor
