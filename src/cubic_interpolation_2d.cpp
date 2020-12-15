#include "corridor/cubic_spline/cubic_interpolation_2d.h"

#include <limits>

#include "corridor/cubic_spline/cubic_spline_coefficients.h"

namespace corridor {

namespace cubic_spline {

using DynamicMatrix = Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic>;
using DynamicColumnVector = Eigen::Matrix<RealType, Eigen::Dynamic, 1>;

// /////////////////////////////////////////////////////////////////////////////
// Derive Spline coefficients from spline data matrix
// /////////////////////////////////////////////////////////////////////////////

SplineCoefficients2d SplineCoefficientsFromDataMatrix(
    const DataMatrix<RealType>& data) {
  SplineCoefficients2d coefficients;
  coefficients.clear();
  for (DataIdx i = 0, max_idx = (data.cols() - 1); i < max_idx; i++) {
    coefficients.emplace_back(Coefficients2d(data.col(i), data.col(i + 1)));
  }
  return coefficients;
}

// /////////////////////////////////////////////////////////////////////////////
// Natural Spline functions
// /////////////////////////////////////////////////////////////////////////////

DynamicMatrix NaturalSplineMatrixM(const DataRow<RealType>& arc_lengths) {
  const auto size = arc_lengths.cols();
  DynamicMatrix matrix_M;
  matrix_M.resize(size, size);
  matrix_M.setZero();
  for (DataIdx i = 1, max_idx = size - 1; i < max_idx; i++) {
    const auto delta_arc_length_1 = arc_lengths(i) - arc_lengths(i - 1);
    const auto delta_arc_length_2 = arc_lengths(i + 1) - arc_lengths(i);
    matrix_M(i, i - 1) = delta_arc_length_1;
    matrix_M(i, i) = 2.f * (delta_arc_length_1 + delta_arc_length_2);
    matrix_M(i, i + 1) = delta_arc_length_2;
  }
  matrix_M(0, 0) = 1.0;
  matrix_M(size - 1, size - 1) = 1.0;
  return matrix_M;
}

DynamicMatrix NaturalSplineMatrixL(const DataRow<RealType>& arc_lengths) {
  const auto size = arc_lengths.cols();
  DynamicMatrix matrix_L;
  matrix_L.resize(size, size);
  matrix_L.setZero();
  for (DataIdx i = 1, max_idx = size - 1; i < max_idx; i++) {
    const auto delta_arc_length_1 = arc_lengths(i) - arc_lengths(i - 1);
    const auto delta_arc_length_2 = arc_lengths(i + 1) - arc_lengths(i);
    matrix_L(i, i - 1) = 6.f / delta_arc_length_1;
    matrix_L(i, i) = (-6.f / delta_arc_length_1) + (-6.f / delta_arc_length_2);
    matrix_L(i, i + 1) = 6.f / delta_arc_length_2;
  }
  return matrix_L;
}

// /////////////////////////////////////////////////////////////////////////////
// General utility functions
// /////////////////////////////////////////////////////////////////////////////

void ArcLengthApproximation(DataMatrix<RealType>& data,
                            const SplineCoefficients2d& coefficients,
                            const RealType epsilon = g_epsilon_arc_length) {
  assert(static_cast<DataSize>(coefficients.size()) == (data.cols() - 1));
  RealType accumulated_arc_length = 0.f;
  // TODO: simplify the for loop
  for (std::size_t idx = 0, max_idx = coefficients.size(); idx < max_idx;
       idx++) {
    const auto data_idx = static_cast<DataIdx>(idx);

    const RealType delta_arc_length =
        data(kArcLength, data_idx + 1) - data(kArcLength, data_idx);
    data(kArcLength, data_idx) = accumulated_arc_length;

    // new accumulated arc-length for (idx+1)
    accumulated_arc_length += ApproxArclengthGaussLegendre(
        coefficients[idx], delta_arc_length, epsilon);
  }
  // Last point
  data(kArcLength, data.cols() - 1) = accumulated_arc_length;
}

DataMatrix<RealType> NaturalSplineDataMatrixFromPoints(
    const CartesianPoints2D& cartesian_points, const RealType epsilon) {
  // 1) Check integrity of points
  const auto pts_size = cartesian_points.size();
  assert(pts_size > 1);

  // Create spline data matrix which will be filled below
  DataMatrix<RealType> data;
  data.resize(BasicDataTypes::kSize, static_cast<Eigen::Index>(pts_size));

  // Add points to data matrix
  for (std::size_t i = 0; i < pts_size; i++) {
    data.block<2, 1>(kPoint_x, static_cast<DataIdx>(i)) = cartesian_points[i];
  }

  // Initial arc-length calculation
  RealType accumulated_arc_length = 0.f;
  for (DataIdx idx = 0, max = data.cols() - 1; idx < max; idx++) {
    data(kArcLength, idx) = accumulated_arc_length;
    const DataPoint<RealType>& p1 = data.block<2, 1>(kPoint_x, idx);
    const DataPoint<RealType>& p2 = data.block<2, 1>(kPoint_x, idx + 1);
    accumulated_arc_length += ChordLength(p1, p2);
  }
  data(kArcLength, data.cols() - 1) = accumulated_arc_length;

  RealType delta_arc_length = 0.f;
  // Approximation loop for moments and arc-length
  for (int counter = 0; counter < 10; counter++) {
    DynamicMatrix matrixM = NaturalSplineMatrixM(data.row(kArcLength));
    DynamicMatrix matrixL = NaturalSplineMatrixL(data.row(kArcLength));

    // Define Moments
    Eigen::FullPivLU<DynamicMatrix> full_piv_lu(matrixM);
    data.row(kMoment_x) = full_piv_lu.solve(
        matrixL * static_cast<DynamicColumnVector>(data.row(kPoint_x)));
    data.row(kMoment_y) = full_piv_lu.solve(
        matrixL * static_cast<DynamicColumnVector>(data.row(kPoint_y)));

    // Define Coefficients
    SplineCoefficients2d coefficients = SplineCoefficientsFromDataMatrix(data);

    // Define new arc-length
    delta_arc_length = data(kArcLength, data.cols() - 1);
    ArcLengthApproximation(data, coefficients);
    delta_arc_length -= data(kArcLength, data.cols() - 1);
    if (abs(delta_arc_length) < epsilon) {
      break;
    }
  }
  return data;
}

// /////////////////////////////////////////////////////////////////////////////
// Clamped Spline functions
// /////////////////////////////////////////////////////////////////////////////

DynamicMatrix ClampedSplineMatrixM(const DataRow<RealType>& arc_lengths) {
  const auto size = arc_lengths.cols();
  DynamicMatrix matrix_M;
  matrix_M.resize(size, size);
  matrix_M.setZero();
  for (DataIdx i = 1, max_idx = size - 1; i < max_idx; i++) {
    const auto delta_arc_length_1 = arc_lengths(i) - arc_lengths(i - 1);
    const auto delta_arc_length_2 = arc_lengths(i + 1) - arc_lengths(i);
    matrix_M(i, i - 1) = delta_arc_length_1;
    matrix_M(i, i) = 2.f * (delta_arc_length_1 + delta_arc_length_2);
    matrix_M(i, i + 1) = delta_arc_length_2;
  }

  // Set first and last row specifically for the clamped constraints
  const auto first_delta_l = arc_lengths(1) - arc_lengths(0);
  matrix_M(0, 0) = 2 * first_delta_l;
  matrix_M(0, 1) = first_delta_l;

  const auto last_delta_l = arc_lengths(size - 1) - arc_lengths(size - 2);
  matrix_M(size - 1, size - 1) = 2 * last_delta_l;
  matrix_M(size - 1, size - 2) = last_delta_l;
  return matrix_M;
}

DynamicMatrix ClampedSplineMatrixL(const DataRow<RealType>& arc_lengths) {
  const auto row_size = arc_lengths.cols();
  const auto col_size = row_size + 2;
  DynamicMatrix matrix_L;
  matrix_L.resize(row_size,
                  col_size);  // two additional columns for the tangent
  matrix_L.setZero();
  for (DataIdx i = 1, max_idx = row_size - 1; i < max_idx; i++) {
    const auto delta_arc_length_1 = arc_lengths(i) - arc_lengths(i - 1);
    const auto delta_arc_length_2 = arc_lengths(i + 1) - arc_lengths(i);
    matrix_L(i, i) = 6.f / delta_arc_length_1;
    matrix_L(i, i + 1) =
        (-6.f / delta_arc_length_1) + (-6.f / delta_arc_length_2);
    matrix_L(i, i + 2) = 6.f / delta_arc_length_2;
  }

  // Set first and last row specifically for the clamped constraints
  const auto first_delta_l = arc_lengths(1) - arc_lengths(0);
  matrix_L(0, 0) = -6.0;
  matrix_L(0, 1) = -6.0 / first_delta_l;
  matrix_L(0, 2) = 6.0 / first_delta_l;

  const auto last_delta_l =
      arc_lengths(row_size - 1) - arc_lengths(row_size - 2);
  matrix_L(row_size - 1, col_size - 1) = 6.0;
  matrix_L(row_size - 1, col_size - 2) = -6.0 / last_delta_l;
  matrix_L(row_size - 1, col_size - 3) = 6.0 / last_delta_l;

  return matrix_L;
}

DataMatrix<RealType> ClampedSplineDataMatrixFromPoints(
    const CartesianPoints2D& cartesian_points,
    const CartesianVector2D& first_tangent,
    const CartesianVector2D& last_tangent, const RealType epsilon) {
  // 1) Check integrity of points
  const auto pts_size = cartesian_points.size();
  assert(pts_size > 1);

  // Create spline data matrix which will be filled below
  DataMatrix<RealType> data;
  data.resize(BasicDataTypes::kSize, static_cast<Eigen::Index>(pts_size));

  // Add points to data matrix
  for (std::size_t i = 0; i < pts_size; i++) {
    data.block<2, 1>(kPoint_x, static_cast<DataIdx>(i)) = cartesian_points[i];
  }

  // Initial arc-length calculation
  RealType accumulated_arc_length = 0.f;
  for (DataIdx idx = 0, max = data.cols() - 1; idx < max; idx++) {
    data(kArcLength, idx) = accumulated_arc_length;
    const DataPoint<RealType>& p1 = data.block<2, 1>(kPoint_x, idx);
    const DataPoint<RealType>& p2 = data.block<2, 1>(kPoint_x, idx + 1);
    accumulated_arc_length += ChordLength(p1, p2);
  }
  data(kArcLength, data.cols() - 1) = accumulated_arc_length;

  RealType delta_arc_length = 0.f;
  // Approximation loop for moments and arc-length
  for (int counter = 0; counter < 10; counter++) {
    DynamicMatrix matrixM = ClampedSplineMatrixM(data.row(kArcLength));
    DynamicMatrix matrixL = ClampedSplineMatrixL(data.row(kArcLength));

    // Augmented support point vectors with tangent data
    DynamicColumnVector aug_pts_x(pts_size + 2);
    aug_pts_x << first_tangent.x(),
        static_cast<DynamicColumnVector>(data.row(kPoint_x)), last_tangent.x();
    DynamicColumnVector aug_pts_y(pts_size + 2);
    aug_pts_y << first_tangent.y(),
        static_cast<DynamicColumnVector>(data.row(kPoint_y)), last_tangent.y();

    // Define Moments
    Eigen::FullPivLU<DynamicMatrix> full_piv_lu(matrixM);
    data.row(kMoment_x) = full_piv_lu.solve(matrixL * aug_pts_x);
    data.row(kMoment_y) = full_piv_lu.solve(matrixL * aug_pts_y);

    // Define Coefficients
    SplineCoefficients2d coefficients = SplineCoefficientsFromDataMatrix(data);

    // Define new arc-length
    delta_arc_length = data(kArcLength, data.cols() - 1);
    ArcLengthApproximation(data, coefficients);
    delta_arc_length -= data(kArcLength, data.cols() - 1);
    if (abs(delta_arc_length) < epsilon) {
      break;
    }
  }
  return data;
}

}  // namespace cubic_spline
}  // namespace corridor
