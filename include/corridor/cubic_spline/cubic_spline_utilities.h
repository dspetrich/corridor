#pragma once

#include <Eigen/Core>
#include <chrono>
#include <vector>

#include "corridor/basic_types.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/frenet_types.h"

namespace corridor {
namespace cubic_spline {

/**
 * @brief Evaluates the position at a given arc-length on the data segment
 *
 * @param data_segment: spline segment on the cubic spline
 * @param arc_length: local arc-length on spline segment
 * @return Basic2dPoint<RealType: resulting point on spline segments
 */
CartesianPoint2D EvaluatePosition(const DataSegment<RealType>& data_segment,
                                  const RealType arc_length);

/**
 * @brief Evaluates the tangent vector at a given arc-length on the data segment
 *
 * @param data_segment: spline segment on the cubic spline
 * @param arc_length: local arc-length on spline segment
 * @return Basic2dPoint<RealType
 *>: resulting point on spline segments
 */
CartesianPoint2D EvaluateTangent(const DataSegment<RealType>& data_segment,
                                 const RealType arc_length);

/**
 * @brief Evaluates the normal vector at a given arc-length on the data segment
 *
 * @param data_segment: spline segment on the cubic spline
 * @param arc_length: local arc-length on spline segment
 * @return Basic2dPoint<RealType
 *>: resulting normal vetcor on spline segments
 */
CartesianPoint2D EvaluateNormal(const DataSegment<RealType>& data_segment,
                                const RealType arc_length);

/**
 * @brief Evaluates the curvature at a given arc-length on the data segment
 *
 * @param data_segment: spline segment on the cubic spline
 * @param arc_length: local arc-length on spline segment
 * @return Basic2dPoint<RealType
 *>: resulting curvature vectors on spline segments
 */
CartesianPoint2D EvaluateCurvature(const DataSegment<RealType>& data_segment,
                                   const RealType arc_length);

/**
 * @brief Interpolates the curvature value on the given data segment at
 * arc_length. The sign is determined by the cross product between tangent and
 * curvature vector
 *
 * @param data_segment: spline segment on the cubic spline
 * @param arc_length: local arc-length on spline segment
 * @return RealType signed curvature value. Positive for left turns, negative
 * for right turns
 */
RealType InterpolateSignedCurvatureValue(
    const DataSegment<RealType>& data_segment, const RealType arc_length);

/**
 * @brief Return type for creating a frenet frame tripod with signed curvature
 * [0]: Origin
 * [1]: Tangent
 * [2]: Normal
 * [3]: Signed curvature value
 * [4]: Signed curvature change rate value
 * TODO: redo this!
 */
using FrenetFrameTripod = std::tuple<CartesianPoint2D, CartesianPoint2D,
                                     CartesianPoint2D, RealType, RealType>;
FrenetFrameTripod InterpolateFrenetFrameTripod(
    const DataSegment<RealType>& data_segment, const RealType arc_length);

/**
 * @brief Generates all tangent vectors on the spline nodes
 *
 * @param data
 * @return Eigen::Matrix<RealType
 *, 2, Eigen::Dynamic>
 */
Eigen::Matrix<RealType, 2, Eigen::Dynamic> TangentsOnNodes(
    const DataMatrix<RealType>& data);

/**
 * @brief Project a point to tangent at arc-length
 *
 * @param data_segment
 * @param segment_info
 * @param point
 * @return RealType
 *: arc-length along segment
 */
RealType TangentialProjection(const CartesianPoint2D& point,
                              const Coefficients2d& segment_coeffs,
                              const RealType arc_length);

RealType TangentialProjectionNewtonRaphson(const CartesianPoint2D& point,
                                           const Coefficients2d& segment_coeffs,
                                           const RealType arc_length);

/**
 * @brief FindSegmentCandidates: returns all segment ids which are likely to
 *        contain the projection of the point. In the case that no
 * pependicular projection onto the spline can be found, the point is either
 * befor the first or after the last segment (or both, in case of a circle).
 *
 * @param data: spline data matrix
 * @param point: which is projected on the spline
 * @param segment_points: filled with the found segment candidates
 * @return true: if a perpendicular match is found
 * @return false: if no perpendicular match is found. In this case either
 * the first or last segment or both segments are returned.
 */
bool FindSegmentCandidates(
    const DataMatrix<RealType>& data, const CartesianPoint2D& point,
    SegmentInfoVector<DataIdx, RealType>* segment_points);

/**
 * @brief
 *
 * @param data_segment
 * @param segment_info
 * @param point
 * @param epsilon
 * @return true
 * @return false
 */
bool FindProjectionOnSegment(const DataSegment<RealType>& data_segment,
                             SegmentInfo<DataIdx, RealType>* segment_info,
                             const CartesianPoint2D& point,
                             const RealType epsilon = g_epsilon_projection);

/**
 * @brief
 *
 * @param data_segment
 * @param segment_info
 * @return FrenetFrame2D
 */
FrenetFrame2D ConstructFrenetFrame(const DataSegment<RealType>& data_segment,
                                   SegmentInfo<DataIdx, RealType> segment_info,
                                   const IdType id = InvalidId);

/**
 * @brief
 *
 * @param data
 * @param point
 * @return FrenetFrames2D
 */
FrenetFrames2D ConstructFrenetFrames(const DataMatrix<RealType>& data,
                                     const CartesianPoint2D& point);

/**
 * @brief
 *
 * @param data
 * @param point
 * @param id
 * @return FrenetPositionsWithFrames
 */
FrenetPositionsWithFrames ConstructFrenetPositionsWithFrames(
    const DataMatrix<RealType>& data, const CartesianPoint2D& point,
    const IdType id = InvalidId);

/**
 * @brief Converts list of x,y points to a FrenetPolyline
 *
 * @param data
 * @param x_vec
 * @param y_vec
 * @return FrenetPolyline
 */
FrenetPolyline ConvertToFrenetPolyline(const DataMatrix<RealType>& data,
                                       const CartesianPoints2D& pts);

}  // namespace cubic_spline
}  // namespace corridor
