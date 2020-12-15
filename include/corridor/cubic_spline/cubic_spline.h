#pragma once

#include <Eigen/Core>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <vector>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline_types.h"
#include "corridor/frenet_types.h"

namespace corridor {
namespace cubic_spline {

class CubicSpline {
 public:
  CubicSpline(const IdType id = InvalidId,
              const RealType epsilon = g_epsilon_arc_length)
      : id_(id), epsilon_(epsilon) {}
  CubicSpline(const IdType id, const CartesianPoints2D& points,
              const RealType epsilon = g_epsilon_arc_length)
      : id_(id), epsilon_(epsilon) {
    bool result = constructSplineData(points);
    if (result == false) {
      // If spline data is corrupted, id is changed to invalid
      id_ = InvalidId;
    }
  }
  CubicSpline(const IdType id, const std::vector<RealType>& x_vec,
              const std::vector<RealType>& y_vec,
              const RealType epsilon = g_epsilon_arc_length)
      : id_(id), epsilon_(epsilon) {
    bool result = constructSplineData(x_vec, y_vec);
    if (result == false) {
      // If spline data is corrupted, id is changed to invalid
      id_ = InvalidId;
    }
  }
  CubicSpline(const IdType id, const CartesianPoints2D& cartesian_points,
              const CartesianVector2D& first_tangent,
              const CartesianVector2D& last_tangent,
              const RealType epsilon = g_epsilon_arc_length) {}

  /** @name Simple public get functions */
  ///@{
  inline IdType GetId() const { return id_; }
  inline DataSize GetSize() const { return data_.cols(); }
  inline RealType GetTotalLength() const {
    return (data_.cols() == 0) ? (0.0) : (data_.rightCols<1>()(kArcLength));
  }
  RealType GetCurvatureAt(const RealType arc_length) const;
  ///@}

  /**
   * @brief FrenetFrames constructs all frenet frames for the perpendicular
   *        projection of point onto the spline
   *
   * @param point: Point, which is projected perpendicular onto the spline
   * @return FrenetFrames2D
   */
  FrenetFrames2D FrenetFrames(CartesianPoint2D point) const;

  FrenetPositionWithFrame getFrenetPositionWithFrame(
      CartesianPoint2D point) const;

  FrenetPolyline toFrenetPolyline(const CartesianPoints2D& points) const;

  // Introspection
  friend std::ostream& operator<<(std::ostream& os, const CubicSpline& cs);

 private:
  /**
   * @brief Constructs the spline data matrix from a list of points
   *
   * @param x_vec: vetcor of x-coordinates
   * @param y_vec: vector of y-coordinates
   * @return true - spline fitting was successfull
   * @return false - spline fitting was not possible
   */
  bool constructSplineData(const std::vector<RealType>& x_vec,
                           const std::vector<RealType>& y_vec);

  /**
   * @brief  Constructs the spline data matrix from a polyline
   *
   * @param points
   * @return true
   * @return false
   */
  bool constructSplineData(const CartesianPoints2D& points);

  /**
   * @brief  Constructs the spline data matrix from a polyline with tangent
   * vectors for the first and last point
   *
   * @param points
   * @return true
   * @return false
   */
  bool constructSplineData(const CartesianPoints2D& points,
                           const CartesianVector2D& first_tangent,
                           const CartesianVector2D& last_tangent);

  /**
   * @brief ToString: Debug information
   *
   * @param print_header
   * @return std::string
   */
  std::string ToString(const bool print_header = true) const;

  IdType id_;                  ///< Unique id of the spline
  RealType epsilon_;           ///< Approximation epsilon of spline fitting
  DataMatrix<RealType> data_;  ///< Data matrix of all sample points information
};                             // namespace cubic_spline

inline std::ostream& operator<<(std::ostream& os, const CubicSpline& cs) {
  using namespace std;
  os << "CubicSpline: " << cs.id_ << " \n";
  os << "IDX\tPT_X\tPT_Y\tM_X\t\tM_Y\t\tARCLEN\n";
  for (int idx = 0; idx < cs.data_.cols(); idx++) {
    os << idx << "\t";
    os << left << setw(8) << setfill(' ') << cs.data_(kPoint_x, idx) << "\t";
    os << left << setw(8) << setfill(' ') << cs.data_(kPoint_y, idx) << "\t";
    os << left << setw(8) << setfill(' ') << cs.data_(kMoment_x, idx) << "\t";
    os << left << setw(8) << setfill(' ') << cs.data_(kMoment_y, idx) << "\t";
    os << left << setw(8) << setfill(' ') << cs.data_(kArcLength, idx) << "\n";
  }
  return os;
}

}  // namespace cubic_spline
}  // namespace corridor
