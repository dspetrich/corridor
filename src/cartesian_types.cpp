#include "corridor/cartesian_types.h"

#include "corridor/basic_types.h"
#include "corridor/unscented_transformation/polar_coordinate_transformation.h"
// /////////////////////////////////////////////////////////////////////////////
// Cartesian State 2D
// /////////////////////////////////////////////////////////////////////////////

namespace corridor {
const PolarStatePtr CartesianState2D::getPolarVelocityStatePtr() {
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

}  // namespace corridor