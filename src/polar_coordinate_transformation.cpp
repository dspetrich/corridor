#include "corridor/unscented_transformation/polar_coordinate_transformation.h"

#include "corridor/basic_types.h"
#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

namespace corridor {
namespace unscented_transformation {

Eigen::Vector2d CartesianToPolarTransformation2D(
    const Eigen::Vector2d& cartesian_vector) {
  Eigen::Vector2d polar_vector;
  const auto& x1 = cartesian_vector(0);
  const auto& x2 = cartesian_vector(1);
  polar_vector(0) = std::sqrt(x1 * x1 + x2 * x2);
  polar_vector(1) = std::atan2(x2, x1);
  return polar_vector;
};

Eigen::Vector2d PolarToCartesianTransformation2D(
    const Eigen::Vector2d& polar_vector) {
  Eigen::Vector2d cartesian_vector;
  const auto& radius = polar_vector(0);
  const auto& phi = polar_vector(1);
  cartesian_vector(0) = radius * std::cos(phi);
  cartesian_vector(1) = radius * std::sin(phi);
  return cartesian_vector;
};

void ToPolarCoordinates2D(const Eigen::Vector2d& initial_x,
                          const Eigen::Matrix2d& initial_P,
                          PolarVector2D* resulting_x,
                          PolarCovarianceMatrix2D* resulting_P) {
  // State tranformation: vel_x, vel_y -> abs_vel, theta
  MerweScaledSigmaPoints<2> sigma_pts_generator;
  const auto& sigmas =
      sigma_pts_generator.generateSigmaPoints(initial_x, initial_P);

  // Transformation function
  Eigen::MatrixXd transformed_sigmas(initial_x.rows(), sigmas.cols());
  for (int i = 0; i < sigmas.cols(); i++) {
    transformed_sigmas.col(i) = CartesianToPolarTransformation2D(sigmas.col(i));
  }

  EstimateStateMeanAndCovarianceMatrix(
      transformed_sigmas, sigma_pts_generator.weightsMean(),
      sigma_pts_generator.weightsCovMat(), (*resulting_x), (*resulting_P), 1);
};

void FromPolarCoordinates2D(const PolarVector2D& initial_x,
                            const PolarCovarianceMatrix2D& initial_P,
                            Eigen::Vector2d* resulting_x,
                            Eigen::Vector2d* resulting_P) {
  // State tranformation: vel_x, vel_y -> abs_vel, theta
  MerweScaledSigmaPoints<2> sigma_pts_generator;

  const auto& sigmas =
      sigma_pts_generator.generateSigmaPoints(initial_x, initial_P);

  // Transformation function
  Eigen::MatrixXd transformed_sigmas(initial_x.rows(), sigmas.cols());
  for (int i = 0; i < sigmas.cols(); i++) {
    transformed_sigmas.col(i) = PolarToCartesianTransformation2D(sigmas.col(i));
  }

  EstimateStateMeanAndCovarianceMatrix(
      transformed_sigmas, sigma_pts_generator.weightsMean(),
      sigma_pts_generator.weightsCovMat(), (*resulting_x), (*resulting_P), 1);
};

}  // namespace unscented_transformation
}  // namespace corridor