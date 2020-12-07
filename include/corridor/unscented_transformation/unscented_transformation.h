#pragma once

#include <Eigen/Core>
#include <cmath>
#include <functional>

#include "corridor/basic_types.h"
#include "corridor/unscented_transformation/sigma_points.h"

namespace corridor {
namespace unscented_transformation {

template <typename SigmaPtsMatrix, typename WeightVector, typename StateVector,
          typename StateCovMatrix>
void EstimateStateMeanAndCovarianceMatrix(
    const Eigen::MatrixBase<SigmaPtsMatrix>& transformed_sigma_pts,
    const Eigen::MatrixBase<WeightVector>& weights_mean,
    const Eigen::MatrixBase<WeightVector>& weights_cov_mat,
    Eigen::MatrixBase<StateVector> const& resulting_mean,
    Eigen::MatrixBase<StateCovMatrix> const& resulting_cov_mat,
    const int angular_value_index = -1) {
  // Mutable eigen-based paramters need to be passed as const
  // reference. To be able to change the value, the constness need to be casted
  // away. More details about this "hack"  and why it is needed can be found
  // here: https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html
  Eigen::MatrixBase<StateVector>& mean =
      const_cast<Eigen::MatrixBase<StateVector>&>(resulting_mean);
  Eigen::MatrixBase<StateCovMatrix>& cov_mat =
      const_cast<Eigen::MatrixBase<StateCovMatrix>&>(resulting_cov_mat);

  // Reset results
  mean.fill(0.0);
  cov_mat.fill(0.0);

  // Initialize state vector and covariance matrix
  const int state_dim = transformed_sigma_pts.rows();
  const int n_sigma_pts = transformed_sigma_pts.cols();

  assert(angular_value_index < state_dim);
  const int j = angular_value_index;  // for easy access

  // ///////////////////////////////////////////////////////////////////////////
  // Estimate state mean
  // ///////////////////////////////////////////////////////////////////////////
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

  // ///////////////////////////////////////////////////////////////////////////
  // Estimate state covariance matrix
  // ///////////////////////////////////////////////////////////////////////////
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

}  // namespace unscented_transformation
}  // namespace corridor