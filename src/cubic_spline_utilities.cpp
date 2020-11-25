#include "corridor/cubic_spline/cubic_spline_utilities.h"

#include "corridor/basic_types.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"

namespace corridor {
namespace cubic_spline {

CartesianPoint2D EvaluatePosition(const DataSegment<RealType>& data_segment,
                                  const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));
  return segment_coeffs.evaluatePosition(arc_length);
}

CartesianPoint2D EvaluateTangent(const DataSegment<RealType>& data_segment,
                                 const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));
  return segment_coeffs.evaluateTangent(arc_length).normalized();
}

CartesianPoint2D EvaluateNormal(const DataSegment<RealType>& data_segment,
                                const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));
  return segment_coeffs.evaluateNormal(arc_length).normalized();
}

CartesianPoint2D EvaluateCurvature(const DataSegment<RealType>& data_segment,
                                   const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));
  return segment_coeffs.evaluateCurvature(arc_length);
}

RealType InterpolateSignedCurvatureValue(
    const DataSegment<RealType>& data_segment, const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));

  const CartesianPoint2D tangent =
      segment_coeffs.evaluateTangent(arc_length).normalized();
  const CartesianPoint2D curvature =
      segment_coeffs.evaluateCurvature(arc_length);

  // Sign of curvature by direction of cross product of tangent and curvature
  // vector
  const RealType sign_cur_value =
      tangent.x() * curvature.y() - tangent.y() * curvature.x();
  const RealType signed_curvature_value =
      (sign_cur_value >= 0.0) ? (curvature.norm()) : (-curvature.norm());

  return signed_curvature_value;
}

FrenetFrameTripod InterpolateFrenetFrameTripod(
    const DataSegment<RealType>& data_segment, const RealType arc_length) {
  const Coefficients2d segment_coeffs(data_segment.col(0), data_segment.col(1));

  const CartesianPoint2D origin = segment_coeffs.evaluatePosition(arc_length);
  // Tangent and Normal are normalized for safety reasons
  const CartesianPoint2D tangent =
      segment_coeffs.evaluateTangent(arc_length).normalized();
  const CartesianPoint2D normal =
      segment_coeffs.evaluateNormal(arc_length).normalized();
  const CartesianPoint2D curvature =
      segment_coeffs.evaluateCurvature(arc_length);
  const CartesianPoint2D ccr = segment_coeffs.evaluateCurvatureChangeRate();

  // Sign of curvature by direction of cross product of tangent and curvature
  // vector
  const RealType sign_cur_value =
      tangent.x() * curvature.y() - tangent.y() * curvature.x();
  const RealType curvature_value =
      (sign_cur_value >= 0.0) ? (curvature.norm()) : (-curvature.norm());

  const RealType sign_ccr_value = tangent.x() * ccr.y() - tangent.y() * ccr.x();
  const RealType ccr_value =
      (sign_ccr_value >= 0.0) ? (ccr.norm()) : (-ccr.norm());

  return std::make_tuple(origin, tangent, normal, curvature_value, ccr_value);
}

Eigen::Matrix<RealType, 2, Eigen::Dynamic> TangentsOnNodes(
    const DataMatrix<RealType>& data) {
  // 1) Construct all tangents
  Eigen::Matrix<RealType, 2, Eigen::Dynamic> tangents;
  tangents.resize(2, data.cols());
  for (DataIdx idx = 0, max_idx = data.cols() - 1; idx < max_idx; idx++) {
    const DataSegment<RealType>& data_segment =
        data.block<kSize, 2>(kPoint_x, idx);
    tangents.block<2, 1>(0, idx) = EvaluateTangent(data_segment, 0.f);
  }
  // Last tangent
  const DataSegment<RealType>& data_segment =
      data.block<kSize, 2>(kPoint_x, data.cols() - 2);
  const RealType local_l =
      data_segment.col(1)[kArcLength] - data_segment.col(0)[kArcLength];
  tangents.block<2, 1>(0, tangents.cols() - 1) =
      EvaluateTangent(data_segment, local_l);
  return tangents;
}

RealType ProjectPointToTangent(
    const DataSegment<RealType>& data_segment,
    const SegmentInfo<DataIdx, RealType>& segment_info,
    const CartesianPoint2D& point) {
  const CartesianPoint2D origin =
      EvaluatePosition(data_segment, segment_info.relative_arc_length);
  const CartesianPoint2D tangent =
      EvaluateTangent(data_segment, segment_info.relative_arc_length);
  return tangent.dot(point - origin);
}

bool FindSegmentCandidates(
    const DataMatrix<RealType>& data, const CartesianPoint2D& point,
    SegmentInfoVector<DataIdx, RealType>* segment_points) {
  // 1) all tangents on spline nodes
  Eigen::Matrix<RealType, 2, Eigen::Dynamic> tangents = TangentsOnNodes(data);

  // 2) relative vectors from spline nodes to point
  const auto relative_vectors = -(data.topRows(2).colwise() - point);

  // 3) Project relative vectors onto tangents
  const auto projections = (tangents.transpose() * relative_vectors).diagonal();

  // 4) Collect all segments where the projection on the start node is
  // positive and negative for the end node
  for (DataIdx idx = 0, n = projections.rows() - 1; idx < n; idx++) {
    if (0.0 <= projections[idx] && projections[idx + 1] < 0.f) {
      segment_points->emplace_back(
          SegmentInfo<DataIdx, RealType>(idx, projections[idx]));
    }
  }

  // 5) Check if point is located before first node and/or behind last node
  if (segment_points->empty()) {
    if (projections[0] <= 0.f) {
      // Point before first node
      segment_points->emplace_back(SegmentInfo<DataIdx, RealType>(0, 0.0));
    }

    if (0.0 <= projections[projections.rows() - 1]) {
      // Point behind last node
      const auto last_node_idx = projections.rows() - 1;
      const auto delta_s = (data.col(last_node_idx)[kArcLength] -
                            data.col(last_node_idx - 1)[kArcLength]);
      const auto last_segment_idx = last_node_idx - 1;
      segment_points->emplace_back(
          SegmentInfo<DataIdx, RealType>(last_segment_idx, delta_s));
    }
    return false;
  }
  return true;
}

bool FindProjectionOnSegment(const DataSegment<RealType>& data_segment,
                             SegmentInfo<DataIdx, RealType>* segment_info,
                             const CartesianPoint2D& point,
                             const RealType epsilon) {
  RealType l_min = 0.f;
  RealType l_max =
      data_segment.col(1)[kArcLength] - data_segment.col(0)[kArcLength];
  RealType delta_l = 0.f;

  // Check current arc-length in segment_info
  if (segment_info->relative_arc_length <= l_min) {
    segment_info->relative_arc_length = 0.f;
    return true;
  }
  if (l_max < segment_info->relative_arc_length) {
    segment_info->relative_arc_length = l_max;
    delta_l = ProjectPointToTangent(data_segment, *segment_info, point);
    if (std::abs(delta_l) < epsilon) {
      return true;
    }
    if (0.0 < delta_l) {
      return false;
    }
  }

  //! 1) Binary Search (find the closest base point)
  for (int count = 0; count < 20; count++) {
    segment_info->relative_arc_length = 0.5f * (l_max - l_min) + l_min;
    delta_l = ProjectPointToTangent(data_segment, *segment_info, point);

    if (std::abs(delta_l) < epsilon) {
      return true;
    } else if (delta_l > 0.f) {
      l_min = segment_info->relative_arc_length;
    } else {
      l_max = segment_info->relative_arc_length;
    }
  }

  // Reset l_min and l_max
  l_min = 0.f;
  l_max = data_segment.col(1)[kArcLength] - data_segment.col(0)[kArcLength];

  // Newton's method base projection search
  for (int count = 0; count < 100; count++) {
    segment_info->relative_arc_length =
        limit(segment_info->relative_arc_length + delta_l, l_min, l_max);
    delta_l = ProjectPointToTangent(data_segment, *segment_info, point);

    if (std::abs(delta_l) < epsilon) {
      return true;
    }
  }
  return false;
}

FrenetFrame2D ConstructFrenetFrame(const DataSegment<RealType>& data_segment,
                                   SegmentInfo<DataIdx, RealType> segment_info,
                                   const IdType id) {
  const RealType arc_length =
      data_segment.col(0)[kArcLength] + segment_info.relative_arc_length;
  // Fill Frenet Frame structure
  const SegmentInfo<IdxType, RealType> seg_point(
      static_cast<IdxType>(segment_info.idx), segment_info.relative_arc_length);

  const FrenetFrameTripod tuple = InterpolateFrenetFrameTripod(
      data_segment, segment_info.relative_arc_length);

  const CartesianPoint2D& origin = std::get<0>(tuple);
  const CartesianPoint2D& tangent = std::get<1>(tuple);
  const CartesianPoint2D& normal = std::get<2>(tuple);
  const RealType curvature = std::get<3>(tuple);
  const RealType curvature_change_rate = std::get<4>(tuple);

  const RealType orientation = std::atan2(tangent.y(), tangent.x());

  const FrenetBase2D frenet_base(id, arc_length, orientation, curvature,
                                 curvature_change_rate, seg_point);
  const FrenetFrame2D frenet_frame(frenet_base, origin, tangent, normal);

  return frenet_frame;
}

FrenetFrames2D ConstructFrenetFrames(const DataMatrix<RealType>& data,
                                     const CartesianPoint2D& point) {
  FrenetPositionsWithFrames positions_with_frames =
      ConstructFrenetPositionsWithFrames(data, point);
  FrenetFrames2D frenet_frames;
  std::transform(positions_with_frames.begin(), positions_with_frames.end(),
                 std::back_inserter(frenet_frames),
                 [](const FrenetPositionWithFrame& p) { return p.frame; });
  return frenet_frames;
}

FrenetPositionsWithFrames ConstructFrenetPositionsWithFrames(
    const DataMatrix<RealType>& data, const CartesianPoint2D& point,
    const IdType id) {
  // Get all relevant spline segments (perpendicular projection)
  SegmentInfoVector<DataIdx, RealType> segment_candidates;
  bool matched = FindSegmentCandidates(data, point, &segment_candidates);

  // Find pependicular projection of point onto segments
  FrenetPositionsWithFrames positions_with_frames;
  for (auto& segment : segment_candidates) {
    const DataSegment<RealType>& data_segment =
        data.block<kSize, 2>(kPoint_x, segment.idx);
    if (matched) {
      // Only try to find the projection point if the point was matched to the
      // segments. In case it wasn't matched, the point is located before or
      // after the spline (or both, in case of a circle like spline). Then,
      // the arc-length is already correct.
      const bool success =
          FindProjectionOnSegment(data_segment, &segment, point);
      if (!success) {
        std::cout << "ERROR: no valid Frenet Frame found on segment candidate: "
                  << segment.idx << std::endl;
        assert(false);
      }
    }
    positions_with_frames.emplace_back(
        ConstructFrenetFrame(data_segment, segment, id), point);
  }
  return positions_with_frames;
}

FrenetPoint2D ConvertToFrenetPoint2D(const DataMatrix<RealType>& data,
                                     const CartesianPoint2D& point) {
  // Construct frenet frame candidates for conversion
  FrenetPositionsWithFrames positions_with_frames =
      ConstructFrenetPositionsWithFrames(data, point);
  // Get frenet point with smallest deviation from the reference line
  FrenetPositionsWithFrames::iterator p_iter = std::min_element(
      positions_with_frames.begin(), positions_with_frames.end(),
      [](const FrenetPositionWithFrame& a, const FrenetPositionWithFrame& b) {
        return a.position.d_value() < b.position.d_value();
      });
  assert(p_iter != positions_with_frames.end());
  return p_iter->position;
}

FrenetPolyline ConvertToFrenetPolyline(const DataMatrix<RealType>& data,
                                       const CartesianPoints2D& cartesian_pts) {
  FrenetPolyline polyline(cartesian_pts.size());
  for (int i = 0, n = cartesian_pts.size(); i < n; i++) {
    polyline.SetPoint(i, ConvertToFrenetPoint2D(data, cartesian_pts[i]));
  }
  return polyline;
}

}  // namespace cubic_spline
}  // namespace corridor
