#include "corridor/frenet_types.h"

#include <Eigen/Core>
#include <iostream>

#include "corridor/cubic_spline/cubic_spline_types.h"

using namespace corridor;

// /////////////////////////////////////////////////////////////////////////////
// Frenet Frame
// /////////////////////////////////////////////////////////////////////////////

CartesianPoint2D FrenetFrame2D::FromFrenetPoint(
    const FrenetPoint2D& frenet_point) const {
  // Vector transfromation, translated by the origin of the frenet frame
  return FromFrenetVector(frenet_point) + origin_;
};

CartesianVector2D FrenetFrame2D::FromFrenetVector(
    const FrenetVector2D& frenet_vector) const {
  //! Local point
  const FrenetVector2D relative_vector{
      frenet_vector.l() - frenet_base_.arc_length, frenet_vector.d()};
  // Coordination transformation
  return rotMat_F2C_ * relative_vector;
};

FrenetPoint2D FrenetFrame2D::FromCartesianPoint(
    const CartesianPoint2D& cartesian_position) const {
  //! Define local point
  const CartesianVector2D relative_vector = cartesian_position - origin_;
  FrenetVector2D tmp = FromCartesianVector(relative_vector);
  // Add arc-length of the frenet base
  tmp.l() += frenet_base_.arc_length;
  //! Convert it to a FrenetPoint and return
  return FrenetPoint2D(tmp);
};

FrenetVector2D FrenetFrame2D::FromCartesianVector(
    const CartesianVector2D& cartesian_vector) const {
  //! Coordination transformation
  FrenetVector2D frenet_vector = rotMat_C2F_ * cartesian_vector;
  return frenet_vector;
};

RealType FrenetFrame2D::FromCartesianOrientation(
    const RealType cartesian_orientation) const {
  const auto angle = cartesian_orientation - frenet_base_.orientation;
  return constrainAngle(angle);
}

// /////////////////////////////////////////////////////////////////////////////
// Frenet Polyline
// /////////////////////////////////////////////////////////////////////////////

RealType FrenetPolyline::deviationAt(const RealType query_l) const {
  // Check if query l value is smaller or larger that the polyline
  const auto s_min = data_(DataType::kArclength, 0);
  const auto s_max = data_(DataType::kArclength, data_.cols() - 1);

  if (query_l <= (s_min + cubic_spline::g_epsilon_projection)) {
    return data_(DataType::kDeviation, 0);
  }
  if ((s_max - cubic_spline::g_epsilon_projection) <= query_l) {
    return data_(DataType::kDeviation, data_.cols() - 1);
  }

  //! Find index of closest data point to a given l value.
  auto start_index = 0;
  (data_.colwise() - FrenetPoint2D(query_l, 0))
      .topRows<1>()
      .colwise()
      .lpNorm<1>()
      .minCoeff(&start_index);

  //! If query_l is smaller then the l-value of the closest node, shift the
  //! start index by one
  const auto s_node = data_(DataType::kArclength, start_index);
  if (start_index > 0 && query_l < s_node) {
    start_index--;
  }

  if (start_index + 1 >= data_.cols()) {
    return data_(DataType::kDeviation, start_index);
  }

  const auto delta_l = data_(DataType::kArclength, start_index + 1) -
                       data_(DataType::kArclength, start_index);
  assert(abs(delta_l) > 1e-3);
  const auto alpha =
      (query_l - data_(DataType::kArclength, start_index)) / delta_l;

  const auto delta_d = data_(DataType::kDeviation, start_index + 1) -
                       data_(DataType::kDeviation, start_index);

  return delta_d * alpha + data_(DataType::kDeviation, start_index);
};