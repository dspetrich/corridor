#pragma once

#include <boost/math/distributions/normal.hpp>
#include <cmath>

#include "corridor/basic_types.h"

namespace corridor {

/**
 * @brief Simple oriented bounding box
 *
 * This OBB class is special in the sense, that it doesn't carry a (center)
 * position. For all operations, where a center is needed it has to be passed
 * additionally. This is done to not store the position information twice with
 * the risk of diverging.
 *
 */
class OrientedBoundingBox {
 public:
  OrientedBoundingBox(const RealType orientation = 0.0,
                      const RealType length = 0.0, const RealType width = 0.0)
      : orientation_(orientation), length_(length), width_(width) {}

  OrientedBoundingBox(const UncertainValue &orientation,
                      const UncertainValue &length, const UncertainValue &width)
      : orientation_(orientation), length_(length), width_(width) {}

  OrientedBoundingBox(const OrientedBoundingBox &other)
      : orientation_(other.orientation_),
        length_(other.length_),
        width_(other.width_) {}

  // getter
  const UncertainValue &orientation() const { return orientation_; };
  const UncertainValue &length() const { return length_; };
  const UncertainValue &width() const { return width_; };

  /**
   * @brief Projection of the shape onto the x1 and x2 axis in which the
   * orientation angle is given. It is assumed that the orientation is
   * defined  mathematically with respect to the x1 axis.
   *
   * @return std::pair<RealType, RealType>: first: projection onto x1 axis
   *                                        second: projection onto x2 axis
   */
  std::pair<RealType, RealType> projection() const {
    const auto cos = std::cos(orientation_.value);
    const auto sin = std::sin(orientation_.value);
    Eigen::Matrix2d projectionMatrix;
    projectionMatrix << cos, sin, sin, cos;
    const Eigen::Vector2d dimension(length_.value, width_.value);
    const Eigen::Vector2d projection = projectionMatrix * dimension;
    return std::make_pair(std::abs(projection.x()), std::abs(projection.y()));
  }

  RealType projectionX1() const {
    const Eigen::Vector2d projection(std::cos(orientation_.value),
                                     std::sin(orientation_.value));
    const Eigen::Vector2d dimension(length_.value, width_.value);
    return std::abs(projection.dot(dimension));
  }

  RealType projectionX2() const {
    const Eigen::Vector2d projection(std::sin(orientation_.value),
                                     std::cos(orientation_.value));
    const Eigen::Vector2d dimension(length_.value, width_.value);
    return std::abs(projection.dot(dimension));
  }

  UncertainValue enclosingCircleRadius() const noexcept {
    const RealType squared_length = length_.value * length_.value;
    const RealType squared_width = width_.value * width_.value;
    const RealType squared_value = 0.25 * (squared_length + squared_width);

    // Assuming length_ and width_ are uncorrelated
    const RealType variance =
        0.25 * ((squared_length / squared_value) * length_.variance() +
                (squared_width / squared_value) * width_.variance());

    return {std::sqrt(squared_value), variance};
  }

 private:
  UncertainValue orientation_;  // [-pi,pi]
  UncertainValue length_;
  UncertainValue width_;
};
// introspection
inline std::ostream &operator<<(std::ostream &os,
                                const OrientedBoundingBox &obb) {
  os << "OBB:\n\torientation: " << obb.orientation();
  os << "\tlength: " << obb.length();
  os << "\twidth: " << obb.width() << "\n";
  return os;
};

}  // namespace corridor
