#pragma once

#include <Eigen/Core>
#include <cmath>
#include <functional>

#include "corridor/basic_types.h"
#include "corridor/unscented_transformation/sigma_points.h"

namespace corridor {
namespace unscented_transformation {

void EstimateStateMeanAndCovarianceMatrix(
    const Eigen::MatrixXd& transformed_sigma_pts,
    const Eigen::VectorXd& weights_mean, const Eigen::VectorXd& weights_cov_mat,
    Eigen::VectorXd* resulting_mean, Eigen::MatrixXd* resulting_cov_mat,
    const int angular_value_index = -1) {
  // Reset results
  resulting_mean->fill(0.0);
  resulting_cov_mat->fill(0.0);

  auto& mean = (*resulting_mean);
  auto& cov_mat = (*resulting_cov_mat);

  // Initialize state vector and covariance matrix
  const int state_dim = transformed_sigma_pts.rows();
  const int n_sigma_pts = transformed_sigma_pts.cols();

  assert(angular_value_index < state_dim);
  const int j = angular_value_index;

  // ////////////////////
  // Estimate state mean
  // ////////////////////
  // Standard formula for mean calculation
  for (int i = 0; i < n_sigma_pts; i++) {  // iterate over sigma points
    mean += weights_mean(i) * transformed_sigma_pts.col(i);
  }

  if (0 <= angular_value_index) {
    // One of the state variables is an angle
    RealType sum_sin = 0.0, sum_cos = 0.0;
    for (int i = 0; i < n_sigma_pts; i++) {  // iterate over sigma points
      sum_sin += weights_mean(i) * std::sin(transformed_sigma_pts(j, i));
      sum_cos += weights_mean(i) * std::cos(transformed_sigma_pts(j, i));
    }
    mean(angular_value_index) = std::atan2(sum_sin, sum_cos);
  }

  // /////////////////////////////////
  // Estimate state covariance matrix
  // /////////////////////////////////
  RealType two_pi = 2 * M_PI;
  for (int i = 0; i < n_sigma_pts; i++) {
    // state difference
    Eigen::VectorXd x_diff = transformed_sigma_pts.col(i) - mean;
    if (0 <= angular_value_index) {
      // constrain angle residual between -pi and pi
      x_diff(j) = std::fmod(x_diff(j), two_pi);
      if (M_PI < x_diff(j)) {
        x_diff(j) -= two_pi;
      }
    }
    cov_mat += weights_cov_mat(i) * x_diff * x_diff.transpose();
  }
};

StateMeanAndCovarianceMatrix EstimateStateMeanAndCovarianceMatrix(
    const Eigen::MatrixXd& transformed_sigma_pts,
    const Eigen::VectorXd& weights_mean, const Eigen::VectorXd& weights_cov_mat,
    const int angular_value_index = -1) {
  // Resulting state
  StateMeanAndCovarianceMatrix state(transformed_sigma_pts.rows());
  EstimateStateMeanAndCovarianceMatrix(transformed_sigma_pts, weights_mean,
                                       weights_cov_mat, &state.mean,
                                       &state.covMat, angular_value_index);
  return state;
};

}  // namespace unscented_transformation
}  // namespace corridor