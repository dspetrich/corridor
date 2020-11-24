#pragma once

#include <Eigen/Core>
#include <vector>

#include "corridor/basic_types.h"

namespace corridor {
namespace cubic_spline {

// Global constant that defines the spline projection approximation uncertainty
constexpr double g_epsilon_projection = 1e-3;
constexpr double g_epsilon_arc_length = 1e-5;

/** Data types of each data point (column) in the data matrix */
enum BasicDataTypes : int8_t {
  kPoint_x = 0,  //!< Cartesian position in x direction [m]
  kPoint_y,      //!< Cartesian position in y direction [m]
  kMoment_x,     //!< Inner moment of the spline in x direction [1/m]
  kMoment_y,     //!< Inner moment of the spline in y direction [1/m]
  kArcLength,    //!< Arc-length from start of spline to specific point [m]
  kSize
};

/**
 * Data Matrix which comprises all needed information of a cubic spline
 * Rows:    sample points  1,2,3,...,i,...,M (M = number of points)
 * Columns: one sample point m_i = [x_i, y_i, m_x_i, m_y_i, l_i, tl_i ]
 */
template <typename T>
using DataMatrix = Eigen::Matrix<T, BasicDataTypes::kSize, Eigen::Dynamic>;
template <typename T>
using DataColumn = Eigen::Matrix<T, BasicDataTypes::kSize, 1>;
template <typename T>
using DataSegment = Eigen::Matrix<T, BasicDataTypes::kSize, 2>;
template <typename T>
using DataRow = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template <typename T>
using DataPoint = Eigen::Matrix<T, 2, 1, Eigen::DontAlign>;
template <typename T>
using DataVector = Eigen::Matrix<T, 2, 1, Eigen::DontAlign>;

using DataIdx = Eigen::Index;
using DataIdxPair = std::pair<DataIdx, DataIdx>;
using DataSize = DataIdx;

/**
 * @brief Segment information
 *
 * @tparam IdType
 * @tparam RType
 */
template <typename IdType, typename RType>
struct SegmentInfo {
  // Start and End indices of the spline segment start
  IdType idx;                 ///< Segment index on spline
  RType relative_arc_length;  //!< Relative arc-length along segment of a
                              //!< point on the segment

  SegmentInfo(const IdType idx = -1, const RType arc_length = 0.f)
      : idx(idx), relative_arc_length(arc_length) {}
};
template <typename IdType, typename RType>
using SegmentInfoVector = std::vector<SegmentInfo<IdType, RType>>;

template <typename IdType, typename RType>
inline std::ostream& operator<<(
    std::ostream& os, const SegmentInfo<IdType, RType>& segment_info) {
  os << "Segment point: idx = " << segment_info.idx
     << ", arc_length = " << segment_info.relative_arc_length;
  return os;
};

}  // namespace cubic_spline
}  // namespace corridor
