#pragma once

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/corridor.h"
#include "corridor/frenet_types.h"

namespace corridor {
namespace unscented_transformation {

FrenetState2D ToFrenetState(const Corridor& corridor,
                            const CartesianState2D cartesian_state);

}  // namespace unscented_transformation
}  // namespace corridor