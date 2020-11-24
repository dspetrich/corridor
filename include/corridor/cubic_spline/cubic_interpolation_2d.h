#pragma once

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <vector>

#include "corridor/basic_types.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_types.h"
#include "corridor/frenet_types.h"

namespace corridor {
namespace cubic_spline {

/**
 * @brief Creates the spline coefficients for each spline segment from the
 *        spline data matrix:
 *        s_x(l) = a_x + b_x * l + c_x * l*l + d_x * l*l*l
 *        s_y(l) = a_y + b_y * l + c_y * l*l + d_y * l*l*l
 *
 * @param data
 * @return SplineCoefficients2d
 */
SplineCoefficients2d SplineCoefficientsFromDataMatrix(
    const DataMatrix<RealType>& data);

/**
 * @brief Creates a spline data matrix a polyline, where the end points are
 *        curvature (strain-energy) free.
 *
 * @param cartesian_points
 * @param epsilon: allowed error of the spline length approximation
 * @return DataMatrix<RealType>
 */
DataMatrix<RealType> naturalSplineDataMatrixFromPoints(
    const CartesianPoints2D& cartesian_points,
    const RealType epsilon = g_epsilon_arc_length);

}  // namespace cubic_spline
}  // namespace corridor
