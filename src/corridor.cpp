#include "corridor/corridor.h"

#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"

namespace corridor {

namespace cs = corridor::cubic_spline;

// /////////////////////////////////////////////////////////////////////////////
// Corridor
// /////////////////////////////////////////////////////////////////////////////

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const RealType distance_left_boundary,
                   const RealType distance_right_boundary) {
  reference_line_ = cs::CubicSpline(id, reference_line_pts);
  const auto num_pts = reference_line_.GetSize();
  left_bound_ = FrenetPolyline(num_pts);
  right_bound_ = FrenetPolyline(num_pts);
  for (int i = 0; i < num_pts; i++) {
    const auto arc_length = reference_line_.GetArclengthAtIndex(i);
    left_bound_.SetPoint(i, {arc_length, distance_left_boundary});
    right_bound_.SetPoint(i, {arc_length, -distance_right_boundary});
  }
}

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const CartesianVector2D& first_tangent,
                   const CartesianVector2D& last_tangent,
                   const RealType distance_left_boundary,
                   const RealType distance_right_boundary) {
  reference_line_ =
      cs::CubicSpline(id, reference_line_pts, first_tangent, last_tangent);
  const auto num_pts = reference_line_.GetSize();
  left_bound_ = FrenetPolyline(num_pts);
  right_bound_ = FrenetPolyline(num_pts);
  for (int i = 0; i < num_pts; i++) {
    const auto arc_length = reference_line_.GetArclengthAtIndex(i);
    left_bound_.SetPoint(i, {arc_length, distance_left_boundary});
    right_bound_.SetPoint(i, {arc_length, -distance_right_boundary});
  }
}  // namespace corridor

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const CartesianPoints2D& left_boundary_pts,
                   const CartesianPoints2D& right_boundary_pts) {
  reference_line_ = cs::CubicSpline(id, reference_line_pts);
  left_bound_ = reference_line_.toFrenetPolyline(left_boundary_pts);
  right_bound_ = reference_line_.toFrenetPolyline(right_boundary_pts);
}

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const CartesianVector2D& first_tangent,
                   const CartesianVector2D& last_tangent,
                   const CartesianPoints2D& left_boundary_pts,
                   const CartesianPoints2D& right_boundary_pts) {
  reference_line_ =
      cs::CubicSpline(id, reference_line_pts, first_tangent, last_tangent);
  left_bound_ = reference_line_.toFrenetPolyline(left_boundary_pts);
  right_bound_ = reference_line_.toFrenetPolyline(right_boundary_pts);
}

BoundaryDistances Corridor::signedDistancesAt(
    const RealType arc_length) const noexcept {
  return std::make_pair(left_bound_.deviationAt(arc_length),
                        right_bound_.deviationAt(arc_length));
}

RealType Corridor::widthAt(const RealType arc_length) const noexcept {
  return left_bound_.deviationAt(arc_length) +
         std::abs(right_bound_.deviationAt(arc_length));
}

RealType Corridor::centerOffset(const RealType arc_length) const noexcept {
  const BoundaryDistances distances = signedDistancesAt(arc_length);
  return (distances.first + distances.second) * 0.5;
}

RealType Corridor::lengthReferenceLine() const noexcept {
  return reference_line_.GetTotalLength();
}

FrenetFrame2D Corridor::FrenetFrame(
    const CartesianPoint2D& position) const noexcept {
  return reference_line_.getFrenetPositionWithFrame(position).frame;
}

FrenetPositionWithFrame Corridor::getFrenetPositionWithFrame(
    const CartesianPoint2D& position,
    const RealType arc_length_hint) const noexcept {
  return reference_line_.getFrenetPositionWithFrame(position, arc_length_hint);
}

std::ostream& operator<<(std::ostream& os, const Corridor& corridor) {
  using namespace std;
  os << "Corridor " << corridor.id() << "\n";
  os << corridor.reference_line_ << "\n";
  os << corridor.left_bound_ << "\n";
  os << corridor.right_bound_ << "\n";
  return os;
}

void Corridor::fillCartesianPolylines(
    const RealType delta_l, CartesianPoints2D* reference_line,
    CartesianPoints2D* left_boundary,
    CartesianPoints2D* right_boundary) const noexcept {
  reference_line->clear();
  left_boundary->clear();
  right_boundary->clear();
  RealType query_l = 0.0;
  RealType max_length = reference_line_.GetTotalLength();
  while (query_l <= max_length) {
    const CartesianPoint2D position = reference_line_.GetPositionAt(query_l);
    const CartesianVector2D normal = reference_line_.GetNormalVectorAt(query_l);
    const RealType d_left = left_bound_.deviationAt(query_l);
    const RealType d_right = right_bound_.deviationAt(query_l);

    reference_line->emplace_back(position);
    left_boundary->emplace_back(position + d_left * normal);
    right_boundary->emplace_back(position + d_right * normal);

    query_l += delta_l;
  }

  if (query_l > max_length) {
    // Add last point
    const CartesianPoint2D position = reference_line_.GetPositionAt(max_length);
    const CartesianVector2D normal =
        reference_line_.GetNormalVectorAt(max_length);
    const RealType d_left = left_bound_.deviationAt(max_length);
    const RealType d_right = right_bound_.deviationAt(max_length);

    reference_line->emplace_back(position);
    left_boundary->emplace_back(position + d_left * normal);
    right_boundary->emplace_back(position + d_right * normal);
  }
}

void Corridor::fillCartesianPolylines(
    CartesianPoints2D* reference_line, CartesianPoints2D* left_boundary,
    CartesianPoints2D* right_boundary) const noexcept {
  reference_line_.fillCartesianPolyline(reference_line);

  left_boundary->clear();
  for (int i = 0, size = left_bound_.size(); i < size; i++) {
    const FrenetPoint2D frenet_point = left_bound_[i];
    const CartesianPoint2D position =
        reference_line_.GetPositionAt(frenet_point.l());
    const CartesianVector2D normal =
        reference_line_.GetNormalVectorAt(frenet_point.l());
    left_boundary->emplace_back(position + frenet_point.d() * normal);
  }

  right_boundary->clear();
  for (int i = 0, size = right_bound_.size(); i < size; i++) {
    const FrenetPoint2D frenet_point = right_bound_[i];
    const CartesianPoint2D position =
        reference_line_.GetPositionAt(frenet_point.l());
    const CartesianVector2D normal =
        reference_line_.GetNormalVectorAt(frenet_point.l());
    right_boundary->emplace_back(position + frenet_point.d() * normal);
  }
}

// /////////////////////////////////////////////////////////////////////////////
// CorridorSequence
// /////////////////////////////////////////////////////////////////////////////

BoundaryDistances CorridorSequence::signedDistancesAt(
    const RealType arc_length) const noexcept {
  auto const_corridor_map_iter = this->get(arc_length);
  const auto delta_arc_length = arc_length - const_corridor_map_iter->first;
  return const_corridor_map_iter->second->signedDistancesAt(delta_arc_length);
}

RealType CorridorSequence::widthAt(const RealType arc_length) const noexcept {
  auto const_corridor_map_iter = this->get(arc_length);
  const auto delta_arc_length = arc_length - const_corridor_map_iter->first;
  return const_corridor_map_iter->second->widthAt(delta_arc_length);
}

RealType CorridorSequence::centerOffsetAt(
    const RealType arc_length) const noexcept {
  auto const_corridor_map_iter = this->get(arc_length);
  const auto delta_arc_length = arc_length - const_corridor_map_iter->first;
  return const_corridor_map_iter->second->centerOffset(delta_arc_length);
}

RealType CorridorSequence::curvatureAt(
    const RealType arc_length) const noexcept {
  auto const_corridor_map_iter = this->get(arc_length);
  const auto delta_arc_length = arc_length - const_corridor_map_iter->first;
  return const_corridor_map_iter->second->curvatureAt(delta_arc_length);
}

RealType CorridorSequence::totalLength() const noexcept {
  return corridor_map_.rbegin()->first /*offset to the last corridor*/ +
         corridor_map_.rbegin()->second->lengthReferenceLine();
}

FrenetPositionWithFrame CorridorSequence::getFrenetPositionWithFrame(
    const CartesianPoint2D& position,
    const RealType start_arc_length) const noexcept {
  CorridorSequenceMap::const_iterator iter = this->get(start_arc_length);
  return this->getFrenetPositionWithFrame(position, iter);
};

FrenetPositionWithFrame CorridorSequence::getFrenetPositionWithFrame(
    const CartesianPoint2D& position,
    CorridorSequenceMap::const_iterator iter) const {
  // Create frenet frame with converted position
  FrenetPositionWithFrame position_with_frame =
      iter->second->getFrenetPositionWithFrame(position);

  // If position is before current corridor, evaluate previous corridor in
  // the sequence
  if (position_with_frame.position.l() < 0.0 && iter != corridor_map_.begin()) {
    iter = std::prev(iter);
    position_with_frame = this->getFrenetPositionWithFrame(position, iter);
    return position_with_frame;
  }

  // If position is behind current corridor, evaluate next corridor in the
  // sequence
  if (iter->second->lengthReferenceLine() < position_with_frame.position.l() &&
      std::next(iter) != corridor_map_.end()) {
    iter = std::next(iter);
    position_with_frame = this->getFrenetPositionWithFrame(position, iter);
    return position_with_frame;
  }

  return position_with_frame;
}

std::ostream& operator<<(std::ostream& os, const CorridorPath& path) {
  using namespace std;
  os << "Corridor-Path:";
  for (const auto& c : path) {
    os << " -> " << c->id();
  }
  os << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const CorridorPaths& paths) {
  using namespace std;
  os << "--- Corridor-Paths ---\n";
  for (const auto& path : paths) {
    os << path << "\n";
  }
  return os;
}

}  // namespace corridor