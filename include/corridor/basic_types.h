#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
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

}  // namespace corridor
