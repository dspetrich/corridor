#pragma once

#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline_types.h"

namespace corridor {
namespace cubic_spline {

struct Coefficients2d {
  RealType a_x, b_x, c_x, d_x;
  RealType a_y, b_y, c_y, d_y;
  RealType max_length;  // coefficients are valid from 0 to max length

  Coefficients2d(const DataColumn<RealType>& q1,
                 const DataColumn<RealType>& q2);

  CartesianPoint2D interpolatePosition(const RealType l) const;
  CartesianPoint2D interpolateTangent(const RealType l,
                                      bool normalized = true) const;
  CartesianPoint2D interpolateNormal(const RealType l) const;
  CartesianPoint2D interpolateCurvature(const RealType l) const;
  CartesianPoint2D interpolateCurvatureChangeRate() const;

  RealType interpolateSignedCurvatureValue(const RealType l) const;
  RealType interpolateSignedCurvatureValue(
      const CartesianPoint2D& tangent, const CartesianPoint2D& curvature) const;

  RealType interpolateSignedCCRValue(const CartesianPoint2D& tangent,
                                     const CartesianPoint2D& ccr) const;

  RealType tangentialProjection(const RealType arc_length,
                                const CartesianPoint2D& point) const;

  RealType tangentialProjectionNewtonRaphson(
      const RealType arc_length, const CartesianPoint2D& point) const;
};

inline std::ostream& operator<<(std::ostream& os,
                                const Coefficients2d& coeffs) {
  os << "x\ty\n";
  os << coeffs.a_x << "\t" << coeffs.a_y << "\n";
  os << coeffs.b_x << "\t" << coeffs.b_y << "\n";
  os << coeffs.c_x << "\t" << coeffs.c_y << "\n";
  os << coeffs.d_x << "\t" << coeffs.d_y << "\n";
  os << std::endl;
  return os;
}

using SplineCoefficients2d = std::vector<Coefficients2d>;

inline CartesianPoint2D Coefficients2d::interpolatePosition(
    const RealType l) const {
  const RealType ll = l * l;
  const RealType lll = ll * l;
  CartesianPoint2D ret_val;
  ret_val << a_x + b_x * l + c_x * ll + d_x * lll,
      a_y + b_y * l + c_y * ll + d_y * lll;
  return ret_val;
}

inline CartesianPoint2D Coefficients2d::interpolateTangent(
    const RealType l, bool normalized) const {
  const RealType ll = l * l;
  CartesianPoint2D ret_val;
  ret_val << (b_x + 2 * c_x * l + 3 * d_x * ll),
      (b_y + 2 * c_y * l + 3 * d_y * ll);
  return (normalized) ? CartesianPoint2D(ret_val.normalized()) : (ret_val);
}

inline CartesianPoint2D Coefficients2d::interpolateNormal(
    const RealType l) const {
  const RealType ll = l * l;
  CartesianPoint2D ret_val;
  ret_val << -(b_y + 2 * c_y * l + 3 * d_y * ll),
      (b_x + 2 * c_x * l + 3 * d_x * ll);
  return ret_val.normalized();
}

inline CartesianPoint2D Coefficients2d::interpolateCurvature(
    const RealType l) const {
  CartesianPoint2D ret_val;
  ret_val << c_x * 2.0 + d_x * l * 6.0, c_y * 2.0 + d_y * l * 6.0;
  return ret_val;
}

inline CartesianPoint2D Coefficients2d::interpolateCurvatureChangeRate() const {
  CartesianPoint2D ret_val;
  ret_val << 6.0 * d_x, 6.0 * d_y;
  return ret_val;
}

inline RealType Coefficients2d::interpolateSignedCurvatureValue(
    const RealType l) const {
  const CartesianPoint2D tangent = interpolateTangent(l);
  const CartesianPoint2D curvature = interpolateCurvature(l);
  return interpolateSignedCurvatureValue(tangent, curvature);
}

inline RealType Coefficients2d::interpolateSignedCurvatureValue(
    const CartesianPoint2D& tangent, const CartesianPoint2D& curvature) const {
  // Sign of curvature by direction of cross product of tangent and curvature
  // vector
  const RealType sign_cur_value =
      tangent.x() * curvature.y() - tangent.y() * curvature.x();
  const RealType signed_curvature_value =
      (sign_cur_value >= 0.0) ? (curvature.norm()) : (-curvature.norm());
  return signed_curvature_value;
}

inline RealType Coefficients2d::interpolateSignedCCRValue(
    const CartesianPoint2D& tangent, const CartesianPoint2D& ccr) const {
  // Sign of ccr by direction of cross product of tangent and ccr
  // vector
  const RealType sign_cur_value = tangent.x() * ccr.y() - tangent.y() * ccr.x();
  const RealType signed_ccr_value =
      (sign_cur_value >= 0.0) ? (ccr.norm()) : (-ccr.norm());
  return signed_ccr_value;
}

inline RealType Coefficients2d::tangentialProjection(
    const RealType arc_length, const CartesianPoint2D& point) const {
  const CartesianPoint2D origin = interpolatePosition(arc_length);
  const CartesianPoint2D tangent = interpolateTangent(arc_length);
  return tangent.dot(point - origin);
}

inline RealType Coefficients2d::tangentialProjectionNewtonRaphson(
    const RealType arc_length, const CartesianPoint2D& point) const {
  // interpolate required curve features at given arc_length
  const CartesianPoint2D origin = interpolatePosition(arc_length);
  const CartesianPoint2D tangent = interpolateTangent(arc_length);
  const CartesianPoint2D normal = interpolateNormal(arc_length);
  const CartesianPoint2D curvature = interpolateCurvature(arc_length);
  const RealType signed_curvature =
      interpolateSignedCurvatureValue(tangent, curvature);

  const CartesianPoint2D delta = point - origin;
  const RealType projection = tangent.dot(delta);
  const RealType projection_derivative =
      signed_curvature * normal.dot(delta) - 1.0;

  if (projection_derivative == 0.0) {
    return 0.0;
  }
  return -(projection / projection_derivative);
}

// /////////////////////////////////////////////////////////////////////////////
// Arc-length approximation utilities
// /////////////////////////////////////////////////////////////////////////////

inline RealType ChordLength(const CartesianPoint2D& p1,
                            const CartesianPoint2D& p2) {
  return (p2 - p1).norm();
}

inline RealType ChordLength(const RealType p1x, const RealType p1y,
                            const RealType p2x, const RealType p2y) {
  BasicVector2D tmp_vector;
  tmp_vector << p2x - p1x, p2y - p1y;
  return tmp_vector.norm();
}

RealType ApproxArclengthNewtonCotes(const Coefficients2d& coeffs,
                                    const RealType approx_h);

RealType ApproxArclengthGaussLegendre(const Coefficients2d& coeffs,
                                      const RealType approx_length,
                                      const RealType epsilon);

RealType SubdivideIntegrationGaussLegendre(const Coefficients2d& coeffs,
                                           const RealType left,
                                           const RealType right,
                                           const RealType approximated_interval,
                                           const RealType epsilon);

RealType IntegralLegendreGauss(const Coefficients2d& coeffs,
                               const RealType left, const RealType right);

}  // namespace cubic_spline
}  // namespace corridor
