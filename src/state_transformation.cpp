#include "corridor/unscented_transformation/state_transformation.h"

#include "corridor/basic_types.h"
#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

namespace corridor {
namespace unscented_transformation {

FrenetState2D ToFrenetState(const Corridor& corridor,
                            const CartesianState2D cartesian_state,
                            const bool moving_frenet_frame) {
  // State transformation: r_x, r_y, vel_x, vel_y -> r_l, r_d, vel_l, vel_d
  MerweScaledSigmaPoints<4> sigma_pts_generator(1.0);

  const auto& sigmas = sigma_pts_generator.generateSigmaPoints(
      cartesian_state.mean(), cartesian_state.covarianceMatrix());

  // Transformation function
  MerweScaledSigmaPoints<4>::SigmaPtsMatrixType transformed_sigmas;
  for (int i = 0; i < sigmas.cols(); i++) {
    const auto frenet_frame = corridor.FrenetFrame(sigmas.block<2, 1>(0, i));
    transformed_sigmas.col(i) = frenet_frame.FromCartesianStateVector(
        sigmas.col(i), moving_frenet_frame);
  }

  FrenetState2D frenet_state;
  EstimateStateMeanAndCovarianceMatrix(
      transformed_sigmas, sigma_pts_generator.weightsMean(),
      sigma_pts_generator.weightsCovMat(), frenet_state.mean(),
      frenet_state.covarianceMatrix());

  return frenet_state;
}

}  // namespace unscented_transformation
}  // namespace corridor