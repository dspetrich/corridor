#pragma once

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace corridor {

// /////////////////////////////////////////////////////////////////////////////
// Basic types
// /////////////////////////////////////////////////////////////////////////////

using IdType = int64_t;
constexpr IdType InvalId = -1;
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

}  // namespace corridor
