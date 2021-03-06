#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iomanip>

#include "corridor/basic_types.h"

namespace corridor {
namespace unscented_transformation {

/**
 * @brief Generates sigma points and weights according to Van der Merwe's 2004
 * dissertation for the Unscented Kalman Filter.
 * It parametrizes the sigma points using alpha, beta, kappa terms, and is the
 * version seen in most publications.
 *
 * @tparam NDim: defines the dimension of the state vector, sigma points and
 * weights
 */
template <int NDim>
class MerweScaledSigmaPoints {
  // Constant definitions
  static constexpr int n_ = NDim;
  static constexpr int n_sigma_pts_ = 2 * n_ + 1;

 public:
  using WeightVectorType =
      Eigen::Matrix<RealType, n_sigma_pts_, 1, Eigen::DontAlign>;
  using SigmaPtsMatrixType =
      Eigen::Matrix<RealType, n_, n_sigma_pts_, Eigen::DontAlign>;
  using CovarianceMatrixType =
      Eigen::Matrix<RealType, n_, n_, Eigen::DontAlign>;

  MerweScaledSigmaPoints(const RealType alpha = 1e-3, const RealType beta = 2.0,
                         const RealType kappa = 0.0)
      : alpha_(alpha), beta_(beta), kappa_(kappa) {
    // Determine lambda
    RealType float_n = static_cast<RealType>(n_);
    lambda_ = alpha_ * alpha_ * (float_n + kappa_) - float_n;
    sqrt_lambda_n_ = std::sqrt(lambda_ + float_n);

    initializeWeights();
  }

  template <typename DerivedVector, typename DerivedMatrix>
  const SigmaPtsMatrixType& generateSigmaPoints(
      const Eigen::MatrixBase<DerivedVector>& x,
      const Eigen::MatrixBase<DerivedMatrix>& P) {
    assert(x.rows() == n_ && P.rows() == n_ && P.cols() == n_);

    // create square root matrix
    const CovarianceMatrixType sqrt_P = P.llt().matrixL();

    sigma_pts_.col(0) = x;
    for (int i = 0; i < n_; i++) {
      const auto delta_x = sqrt_lambda_n_ * sqrt_P.col(i);
      // columns 1 -> n_ = x + sqrt((lambda + n_aug_) * P_)
      sigma_pts_.col(i + 1) = x + delta_x;
      // columns n_aug_+1 -> 2*n_ = x - sqrt((lambda + n_aug_) * P_)
      sigma_pts_.col(i + 1 + n_) = x - delta_x;
    }

    return sigma_pts_;
  }

  const WeightVectorType& weightsMean() const { return weights_mean_; };
  const WeightVectorType& weightsCovMat() const { return weights_cov_; };

 private:
  // Determins the spread of the sigma points around the mean.
  RealType alpha_;

  // beta incorporates prior knowledge of the distribution of the mean. For
  // Gaussian x beta=2 is optimal
  RealType beta_;

  // Secondary scaling parameter usually set to 0 or (3-n)
  RealType kappa_;

  // Resulting sigma point spreading parameter
  RealType lambda_;

  RealType sqrt_lambda_n_;

  WeightVectorType weights_mean_;
  // Weight for each sigma point for the covariance matrix calculation
  WeightVectorType weights_cov_;

  SigmaPtsMatrixType sigma_pts_;

  void initializeWeights() {
    const RealType c = 0.5 / (RealType(n_) + lambda_);
    weights_mean_.fill(c);
    weights_cov_.fill(c);

    // weight for covariance matrix's first point is different than for the
    // mean
    weights_mean_[0] = (lambda_) / (lambda_ + n_);
    weights_cov_[0] = weights_mean_[0] + (1 - (alpha_ * alpha_) + beta_);
  }
};

}  // namespace unscented_transformation
}  // namespace corridor