#include "corridor/cubic_spline/cubic_spline_coefficients.h"

namespace corridor {
namespace cubic_spline {

Coefficients2d::Coefficients2d(const DataColumn<RealType>& q1,
                               const DataColumn<RealType>& q2) {
  const DataColumn<RealType> delta_q = q2 - q1;
  a_x = q1(kPoint_x);
  a_y = q1(kPoint_y);

  b_x = delta_q(kPoint_x) / delta_q(kArcLength) -
        (delta_q(kArcLength) / 6.f) * (2 * q1(kMoment_x) + q2(kMoment_x));

  b_y = delta_q(kPoint_y) / delta_q(kArcLength) -
        (delta_q(kArcLength) / 6.f) * (2 * q1(kMoment_y) + q2(kMoment_y));
  c_x = 0.5 * q1(kMoment_x);
  c_y = 0.5 * q1(kMoment_y);

  d_x = delta_q(kMoment_x) / (6.0 * delta_q(kArcLength));
  d_y = delta_q(kMoment_y) / (6.0 * delta_q(kArcLength));
}

RealType ApproxArclengthNewtonCotes(const Coefficients2d& coeffs,
                                    const RealType approx_h) {
  static RealType weight[7] = {41.f, 216.f, 27.f, 272.f, 27.f, 216.f, 41.f};

  RealType new_h = 0.0;
  RealType u = 0.0;

  BasicVector2D t_vec;
  for (int k = 0; k < 7; k++) {
    u = static_cast<RealType>(k / 6.f * approx_h);

    t_vec = coeffs.evaluateTangent(u);

    new_h += weight[k] / 840.f * t_vec.norm();
  }

  return new_h * approx_h;
}

RealType ApproxArclengthGaussLegendre(const Coefficients2d& coeffs,
                                      const RealType approx_length,
                                      const RealType epsilon) {
  // define interval of the integral
  RealType left = 0.0;
  RealType right = approx_length;

  // whole integral
  RealType full_interval = IntegralLegendreGauss(coeffs, left, right);

  full_interval = SubdivideIntegrationGaussLegendre(coeffs, left, right,
                                                    full_interval, epsilon);
  return full_interval;
}

RealType SubdivideIntegrationGaussLegendre(const Coefficients2d& coeffs,
                                           const RealType left,
                                           const RealType right,
                                           const RealType approximated_interval,
                                           const RealType epsilon) {
  //! divide interval into two sub halves
  RealType mid = (left + right) * 0.5;

  //! calculate the integral values for both sides
  RealType left_interval = IntegralLegendreGauss(coeffs, left, mid);
  RealType right_interval = IntegralLegendreGauss(coeffs, mid, right);

  RealType full_interval = left_interval + right_interval;

  // if accuracy is not good enough, divide the integration interval
  if (fabs(approximated_interval - full_interval) > epsilon) {
    RealType left_sub = SubdivideIntegrationGaussLegendre(
        coeffs, left, mid, left_interval, 0.5 * epsilon);
    RealType right_sub = SubdivideIntegrationGaussLegendre(
        coeffs, mid, right, right_interval, 0.5 * epsilon);
    return left_sub + right_sub;
  } else {
    return full_interval;
  }
}

RealType IntegralLegendreGauss(const Coefficients2d& coeffs,
                               const RealType left, const RealType right) {
  //! weights
  static const RealType w_2_3 = 5.0 / 9.0;
  static const RealType weights[3] = {8.0 / 9.0, w_2_3, w_2_3};

  //! certain evaluation points
  static const RealType p_2_3 = sqrt(3.0 / 5.0);
  static const RealType p[3] = {0.0, p_2_3, -p_2_3};

  //! integration boundaries
  const RealType var_1 = 0.5 * (right - left);
  const RealType var_2 = 0.5 * (left + right);

  //! result
  RealType res_val = 0.0;
  //! function values
  for (int i = 0; i < 3; i++) {
    res_val += weights[i] * coeffs.evaluateTangent(var_1 * p[i] + var_2).norm();
  }

  return res_val * var_1;
}

}  // namespace cubic_spline
}  // namespace corridor