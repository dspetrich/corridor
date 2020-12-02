#pragma once

#include <Eigen/Core>
#include <cmath>
#include <functional>

#include "corridor/unscented_transformation/sigma_points.h"

namespace corridor {
namespace unscented_transformation {

struct StateMeanAndCovarianceMatrix {
  StateMeanAndCovarianceMatrix(const int state_dim = 0)
      : mean(Eigen::VectorXd(state_dim)),
        covMat(Eigen::MatrixXd(state_dim, state_dim)) {
    mean.fill(0.0);
    covMat.fill(0.0);
  }
  StateMeanAndCovarianceMatrix(const Eigen::VectorXd& _mean,
                               const Eigen::MatrixXd& _covMat)
      : mean(_mean), covMat(_covMat) {}
  Eigen::VectorXd mean;
  Eigen::MatrixXd covMat;
};

// Introspection
std::ostream& operator<<(std::ostream& os,
                         const StateMeanAndCovarianceMatrix& state) {
  os << "StateMeanAndCovarianceMatrix\n";
  os << "mean: " << state.mean.transpose() << "\n";
  os << state.covMat;
  return os;
}

StateMeanAndCovarianceMatrix EstimateStateMeanAndCovarianceMatrix(
    const Eigen::MatrixXd transformed_sigma_pts,
    const Eigen::VectorXd weights_mean, const Eigen::VectorXd weights_cov_mat,
    const int angular_value_index = -1) {
  // Initialize state vector and covariance matrix
  const int state_dim = transformed_sigma_pts.rows();
  const int n_sigma_pts = transformed_sigma_pts.cols();

  assert(angular_value_index < state_dim);
  const int j = angular_value_index;

  StateMeanAndCovarianceMatrix state(state_dim);

  // ////////////////////
  // Estimate state mean
  // ////////////////////
  // Standard formula for mean calculation
  for (int i = 0; i < n_sigma_pts; i++) {  // iterate over sigma points
    state.mean += weights_mean(i) * transformed_sigma_pts.col(i);
  }

  if (0 <= angular_value_index) {
    // One of the state variables is an angle
    RealType sum_sin = 0.0, sum_cos = 0.0;
    for (int i = 0; i < n_sigma_pts; i++) {  // iterate over sigma points
      sum_sin += weights_mean(i) * std::sin(transformed_sigma_pts(j, i));
      sum_cos += weights_mean(i) * std::cos(transformed_sigma_pts(j, i));
    }
    state.mean(angular_value_index) = std::atan2(sum_sin, sum_cos);
  }

  // /////////////////////////////////
  // Estimate state covariance matrix
  // /////////////////////////////////
  RealType two_pi = 2 * M_PI;
  for (int i = 0; i < n_sigma_pts; i++) {
    // state difference
    Eigen::VectorXd x_diff = transformed_sigma_pts.col(i) - state.mean;
    if (0 <= angular_value_index) {
      // constrain angle residual between -pi and pi
      x_diff(j) = std::fmod(x_diff(j), two_pi);
      if (M_PI < x_diff(j)) {
        x_diff(j) -= two_pi;
      }
    }
    state.covMat += weights_cov_mat(i) * x_diff * x_diff.transpose();
  }

  return state;
};

}  // namespace unscented_transformation
}  // namespace corridor