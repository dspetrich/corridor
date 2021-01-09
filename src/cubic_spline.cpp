#include "corridor/cubic_spline/cubic_spline.h"

#include <iostream>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"

using namespace corridor;
using namespace cubic_spline;

bool CubicSpline::constructSplineData(const std::vector<RealType>& x_vec,
                                      const std::vector<RealType>& y_vec) {
  if (x_vec.size() < 2 || x_vec.size() != y_vec.size()) {
    return false;
  }
  CartesianPoints2D points;
  for (int i = 0, n = x_vec.size(); i < n; i++) {
    points.emplace_back(x_vec[i], y_vec[i]);
  }

  // Natural spline
  data_ = NaturalSplineDataMatrixFromPoints(points, epsilon_);
  return true;
}

bool CubicSpline::constructSplineData(const CartesianPoints2D& points) {
  if (points.size() < 2) {
    return false;
  }
  // Natural spline
  data_ = NaturalSplineDataMatrixFromPoints(points, epsilon_);
  return true;
}

bool CubicSpline::constructSplineData(const CartesianPoints2D& points,
                                      const CartesianVector2D& first_tangent,
                                      const CartesianVector2D& last_tangent) {
  if (points.size() < 2) {
    return false;
  }
  // clamped spline
  data_ = ClampedSplineDataMatrixFromPoints(points, first_tangent, last_tangent,
                                            epsilon_);
  return true;

  return true;
}

CartesianPoint2D CubicSpline::GetPositionAt(const RealType arc_length) const {
  // Get segment index
  DataMatrix<RealType>::Index index = GetSegmentIndexAtArcLength(arc_length);

  // Construct segment data
  const DataSegment<RealType>& data_segment =
      data_.block<kSize, 2>(kPoint_x, index);
  const RealType relative_arc_length = arc_length - data_segment(kArcLength, 0);
  return EvaluatePosition(data_segment, relative_arc_length);
}

CartesianPoint2D CubicSpline::GetNormalVectorAt(
    const RealType arc_length) const {
  // Get segment index
  DataMatrix<RealType>::Index index = GetSegmentIndexAtArcLength(arc_length);

  // Construct segment data
  const DataSegment<RealType>& data_segment =
      data_.block<kSize, 2>(kPoint_x, index);
  const RealType relative_arc_length = arc_length - data_segment(kArcLength, 0);
  return EvaluateNormal(data_segment, relative_arc_length);
}

RealType CubicSpline::GetCurvatureAt(const RealType arc_length) const {
  // Get segment index
  DataMatrix<RealType>::Index index = GetSegmentIndexAtArcLength(arc_length);

  // Construct segment data
  const DataSegment<RealType>& data_segment =
      data_.block<kSize, 2>(kPoint_x, index);

  const RealType relative_arc_length = arc_length - data_segment(kArcLength, 0);

  return InterpolateSignedCurvatureValue(data_segment, relative_arc_length);
}

FrenetFrames2D CubicSpline::FrenetFrames(CartesianPoint2D point) const {
  return ConstructFrenetFrames(data_, point);
}

FrenetPositionWithFrame CubicSpline::getFrenetPositionWithFrame(
    CartesianPoint2D point) const {
  FrenetPositionsWithFrames positions =
      ConstructFrenetPositionsWithFrames(data_, point, id_);

  // Get frenet point with smallest displacement from centerline
  FrenetPositionsWithFrames::iterator p_iter = std::min_element(
      positions.begin(), positions.end(),
      [](const FrenetPositionWithFrame& a, const FrenetPositionWithFrame& b) {
        return a.position.d_value() < b.position.d_value();
      });

  return (*p_iter);
}

FrenetPolyline CubicSpline::toFrenetPolyline(
    const CartesianPoints2D& points) const {
  return ConvertToFrenetPolyline(data_, points);
}

const DataMatrix<RealType>::Index CubicSpline::GetSegmentIndexAtArcLength(
    const RealType arc_length) const noexcept {
  DataMatrix<RealType>::Index index = 0;
  DataMatrix<RealType>::Index max_index = data_.cols() - 2;
  // Get index which is larger than the query arc_length
  const bool valid =
      (data_.row(kArcLength).array() > arc_length).maxCoeff(&index);
  if (!valid) {
    // Arc-length is longer as the refernce line. Return last segment
    return max_index;
  }
  // Check that always a valid data segment is used
  index = (0 < index) ? (index - 1) : (index);
  index = (max_index < index) ? (max_index) : (index);
  return index;
}

std::string CubicSpline::ToString(const bool print_header) const {
  std::ostringstream out;
  out << "CubicSpline: " << id_ << " \n";
  if (print_header) {
    out << "IDX\tPT_X\tPT_Y\tM_X\tM_Y\tARCLEN\n";
  }
  for (int idx = 0; idx < data_.cols(); idx++) {
    out << idx << "\t" << data_(kPoint_x, idx) << "\t" << data_(kPoint_y, idx)
        << "\t" << data_(kMoment_x, idx) << "\t" << data_(kMoment_y, idx)
        << "\t" << data_(kArcLength, idx) << std::endl;
  }

  return out.str();
}