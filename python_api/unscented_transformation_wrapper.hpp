#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

namespace corridor {
namespace unscented_transformation {

struct FlatCartesianStateAndCovMat2D {
  RealType x;
  RealType y;
  RealType var_x;
  RealType var_y;
  RealType cov_xy;
};

struct FlatPolarStateAndCovMat2D {
  RealType r;
  RealType phi;
  RealType var_r;
  RealType var_phi;
  RealType cov_rphi;
};

StateMeanAndCovarianceMatrix Convert(
    const FlatCartesianStateAndCovMat2D& flat_cartesian_state) {
  StateMeanAndCovarianceMatrix state(2);
  state.mean << flat_cartesian_state.x, flat_cartesian_state.y;
  state.covMat(0, 0) = flat_cartesian_state.var_x;
  state.covMat(1, 1) = flat_cartesian_state.var_y;
  state.covMat(1, 0) = flat_cartesian_state.cov_xy;
  state.covMat(0, 1) = flat_cartesian_state.cov_xy;
  return state;
}

FlatPolarStateAndCovMat2D Convert(
    const StateMeanAndCovarianceMatrix& polar_state) {
  FlatPolarStateAndCovMat2D flat_polar_state;
  flat_polar_state.r = polar_state.mean(0);
  flat_polar_state.phi = polar_state.mean(1);
  flat_polar_state.var_r = polar_state.covMat(0, 0);
  flat_polar_state.var_phi = polar_state.covMat(1, 1);
  flat_polar_state.cov_rphi = polar_state.covMat(1, 0);
  return flat_polar_state;
}

Eigen::Vector2d PolarToCartesianTransformation2D(
    const Eigen::Vector2d& polar_vector) {
  Eigen::Vector2d cartesian_vector;
  const auto& radius = polar_vector(0);
  const auto& phi = polar_vector(1);
  cartesian_vector(0) = radius * std::cos(phi);
  cartesian_vector(1) = radius * std::sin(phi);
  return cartesian_vector;
}

Eigen::Vector2d CartesianToPolarTransformation2D(
    const Eigen::Vector2d& cartesian_vector) {
  Eigen::Vector2d polar_vector;
  const auto& x1 = cartesian_vector(0);
  const auto& x2 = cartesian_vector(1);
  polar_vector(0) = std::sqrt(x1 * x1 + x2 * x2);
  polar_vector(1) = std::atan2(x2, x1);
  return polar_vector;
}

StateMeanAndCovarianceMatrix UnscentedTransformationPolarCoordinates(
    const StateMeanAndCovarianceMatrix& initial_state) {
  // State tranformation: vel_x, vel_y -> abs_vel, theta
  MerweScaledSigmaPoints sigma_pts_generator(2);

  const auto& sigmas = sigma_pts_generator.generateSigmaPoints(
      initial_state.mean, initial_state.covMat);

  // Tranformation function
  Eigen::MatrixXd transformed_sigmas(initial_state.mean.size(), sigmas.cols());
  for (int i = 0; i < sigmas.cols(); i++) {
    transformed_sigmas.col(i) = CartesianToPolarTransformation2D(sigmas.col(i));
  }

  return EstimateStateMeanAndCovarianceMatrix(
      transformed_sigmas, sigma_pts_generator.weightsMean(),
      sigma_pts_generator.weightsCovMat(), 1);
}

}  // namespace unscented_transformation
}  // namespace corridor