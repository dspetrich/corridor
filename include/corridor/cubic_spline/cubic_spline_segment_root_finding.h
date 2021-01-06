#pragma once

#include <utility>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_types.h"

namespace corridor {
namespace cubic_spline {

// /////////////////////////////////////////////////////////////////////////////
// Check boundary conditions
// /////////////////////////////////////////////////////////////////////////////
std::tuple<bool, bool, RealType> LimitArcLengthToSegmentLimits(
    const Coefficients2d& coeffs, const RealType arc_length,
    const CartesianPoint2D& point, const RealType epsilon = 1e-3) {
  // Limited arc-length
  RealType limited_arc_length = arc_length;

  // Lower segment bound
  if (arc_length <= 0.0) {
    limited_arc_length = 0.0;
    const RealType delta_l =
        coeffs.tangentialProjection(limited_arc_length, point);
    if (std::abs(delta_l) < epsilon) {
      return std::make_tuple(true, true, limited_arc_length);
    }
    if (delta_l < limited_arc_length) {
      return std::make_tuple(true, false, limited_arc_length);
    }
  }

  // Upper segment bound
  if (coeffs.max_length <= arc_length) {
    limited_arc_length = coeffs.max_length;
    const RealType delta_l =
        coeffs.tangentialProjection(coeffs.max_length, point);
    if (std::abs(delta_l) < epsilon) {
      return std::make_tuple(true, true, coeffs.max_length);
    }
    if (0.0 < delta_l) {
      return std::make_tuple(true, false, coeffs.max_length);
    }
  }
  return std::make_tuple(false, true, limited_arc_length);
}

// /////////////////////////////////////////////////////////////////////////////
// Root-finding algorithms
// /////////////////////////////////////////////////////////////////////////////

/**
 * @brief
 *
 * @param coeffs
 * @param arc_length
 * @param point
 * @param epsilon
 * @return std::pair<bool, RealType>
 */
std::pair<bool, RealType> BisectionMethod(const Coefficients2d& coeffs,
                                          const RealType arc_length,
                                          const CartesianPoint2D& point,
                                          const RealType epsilon = 1e-3) {
  RealType l_min = 0.;
  RealType l_max = coeffs.max_length;
  RealType l_curr = 0.;

  for (unsigned int count = 0; count < 100; count++) {
    const RealType l_curr = 0.5 * (l_max - l_min) + l_min;
    const RealType delta_l = coeffs.tangentialProjection(l_curr, point);

    if (std::abs(delta_l) < epsilon) {
      // std::cout << __FUNCTION__ << ": count " << count
      //           << "; root is: " << l_curr << std::endl;
      return std::make_pair(true, l_curr);
    } else if (delta_l > 0.f) {
      l_min = l_curr;
    } else {
      l_max = l_curr;
    }
  }
  return std::make_pair(false, l_curr);
}

std::pair<bool, RealType> NewtonRaphsonMethod(const Coefficients2d& coeffs,
                                              const RealType arc_length,
                                              const CartesianPoint2D& point,
                                              const RealType epsilon = 1e-3) {
  const RealType l_min = 0.;
  const RealType l_max = coeffs.max_length;
  RealType l_curr = arc_length;
  RealType delta_l = 0.0;

  for (int count = 0; count < 100; count++) {
    l_curr = limit(arc_length + delta_l, l_min, l_max);
    delta_l = coeffs.tangentialProjectionNewtonRaphson(l_curr, point);

    if (std::abs(delta_l) < epsilon) {
      return std::make_pair(true, l_curr);
    }
  }
  return std::make_pair(false, l_curr);
}

/**
 * @brief
 *
 * @param coeffs
 * @param arc_length
 * @param point
 * @param epsilon
 * @return std::pair<bool, RealType>
 */
std::pair<bool, RealType> BrentsMethod(const Coefficients2d& coeffs,
                                       const RealType arc_length,
                                       const CartesianPoint2D& point,
                                       const RealType epsilon = 1e-3) {
  RealType l_min = 0.0;
  RealType l_max = coeffs.max_length;

  RealType f_min = coeffs.tangentialProjection(l_min, point);
  RealType f_max = coeffs.tangentialProjection(l_max, point);

  if (!(f_min * f_max < 0)) {
    // Signs of f_min and f_max must be opposites, root not bracketed
    return std::make_pair(false, 0.0);
  }

  if (std::abs(f_min) < std::abs(f_max)) {
    // if magnitude of f_min is less  than magnitude of f_max
    std::swap(l_min, l_max);
    std::swap(f_min, f_max);
  }

  // Some auxilary variables
  RealType l_curr = l_min;  //<! l_curr now equals the largest magnitude of
                            // the lower and upper bounds
  RealType f_curr = f_min;  //<! precompute function evaluation for point l_curr
                            // by assigning it the same value as f_min

  bool mflag = true;    // boolean flag used to evaluate if statement later on
  RealType l_res = 0.;  // Our Root that will be returned
  RealType f_res = 0.;  // Function value at l_res
  RealType l_d = 0.;    // Only used if mflag is unset (mflag == false)

  for (unsigned int count = 1; count < 100; ++count) {
    // stop if converged on root or error is less than tolerance
    if (std::abs(l_max - l_min) < epsilon) {
      // std::cout << __FUNCTION__ << ": count " << count << "; root is: " <<
      // l_res
      //           << std::endl;
      return std::make_pair(true, l_res);
    }  // end if

    if ((f_min != f_curr) && (f_max != f_curr)) {
      // use inverse quadratic interpolation
      l_res = (l_min * f_max * f_curr / ((f_min - f_max) * (f_min - f_curr))) +
              (l_max * f_min * f_curr / ((f_max - f_min) * (f_max - f_curr))) +
              (l_curr * f_min * f_max / ((f_curr - f_min) * (f_curr - f_max)));
    } else {
      // secant method
      l_res = l_max - f_max * (l_max - l_min) / (f_max - f_min);
    }

    // Crazy BrentsMethod conditions
    // ////////////////////////
    // (condition 1) s is not between  (3f_min+f_max)/4  and f_max
    // (condition 2) (mflag is true and |l_res−l_max| ≥ |l_max−l_curr|/2)
    // (condition 3)(mflag is false and | l_res−l_max | ≥ | l_curr−l_d | / 2)
    // (condition 4) (mflag is set and |l_max−l_curr| < |epsilon|) or
    // (condition 5) (mflag is false and |l_curr−l_d| < |epsilon|)
    if (((l_res < 0.25 * (3 * l_min + l_max)) || (l_max < l_res)) ||
        (mflag &&
         std::abs(l_res - l_max) >= (0.5 * std::abs(l_max - l_curr))) ||
        (!mflag && std::abs(l_res - l_max) >= (0.5 * std::abs(l_curr - l_d))) ||
        (mflag && (std::abs(l_max - l_curr) < epsilon)) ||
        (!mflag && (std::abs(l_curr - l_d) < epsilon))) {
      // bisection method
      l_res = (l_min + l_max) * 0.5;
      mflag = true;
    } else {
      mflag = false;
    }

    f_res = coeffs.tangentialProjection(l_res, point);
    l_d = l_curr;  // first time l_d is set properly

    l_curr = l_max;  // l_curr equal to upper bound
    f_curr = f_max;

    if (f_min * f_res < 0) {
      // f_min and f_res have opposite signs
      l_max = l_res;
      f_max = f_res;
    } else {
      l_min = l_res;
      f_min = f_res;
    }

    if (std::abs(f_min) < std::abs(f_max)) {
      // if magnitude of fa is less than magnitude of fb
      std::swap(l_min, l_max);
      std::swap(f_min, f_max);
    }
  }

  return std::make_pair(false, l_res);
}
}  // namespace cubic_spline
}  // namespace corridor