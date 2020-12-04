#pragma once

#include <boost/math/distributions/normal.hpp>
#include <cmath>

#include "corridor/basic_types.h"

namespace corridor {
namespace math {

RealType evaluateIntegralLineWidthGaussian(const RealType m, const RealType b,
                                           const RealType x,
                                           const RealType sigma_original,
                                           const RealType lower_bound,
                                           const RealType upper_bound) {
  const RealType sigma = (sigma_original <= 1e-30) ? (1e-30) : (sigma_original);

  RealType tau_1 = (lower_bound - x) / sigma;
  RealType tau_2 = (upper_bound - x) / sigma;

  if (std::isnan(tau_1) || std::isnan(tau_2)) {
    return 0.0;
  }

  boost::math::normal n_dist;
  RealType pdf01 = boost::math::pdf(n_dist, tau_1);
  RealType pdf02 = boost::math::pdf(n_dist, tau_2);

  RealType cdf01 = boost::math::cdf(n_dist, tau_1);
  RealType cdf02 = boost::math::cdf(n_dist, tau_2);

  return (m * x + b) * (cdf02 - cdf01) - m * sigma * (pdf02 - pdf01);
}

}  // namespace math
}  // namespace corridor
