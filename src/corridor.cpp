#include "corridor/corridor.h"

namespace corridor {

namespace cs = corridor::cubic_spline;

// /////////////////////////////////////////////////////////////////////////////
// Corridor
// /////////////////////////////////////////////////////////////////////////////

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const RealType distance_left_boundary,
                   const RealType distance_right_boundary) {
  referenceLine_ = cs::CubicSpline(id, reference_line_pts);
  const auto num_pts = referenceLine_.GetSize();
  leftBound_ = FrenetPolyline(num_pts);
  rightBound_ = FrenetPolyline(num_pts);
  for (int i = 0; i < num_pts; i++) {
    const auto arc_length = referenceLine_.GetArclengthAt(i);
    leftBound_.SetPoint(i, {arc_length, distance_left_boundary});
    rightBound_.SetPoint(i, {arc_length, -distance_right_boundary});
  }
}

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const CartesianPoints2D& left_boundary_pts,
                   const CartesianPoints2D& right_boundary_pts) {
  referenceLine_ = cs::CubicSpline(id, reference_line_pts);
  leftBound_ = referenceLine_.toFrenetPolyline(left_boundary_pts);
  rightBound_ = referenceLine_.toFrenetPolyline(right_boundary_pts);
}

Corridor::Corridor(const IdType id, const CartesianPoints2D& reference_line_pts,
                   const CartesianVector2D& first_tangent,
                   const CartesianVector2D& last_tangent,
                   const CartesianPoints2D& left_boundary_pts,
                   const CartesianPoints2D& right_boundary_pts) {
  referenceLine_ =
      cs::CubicSpline(id, reference_line_pts, first_tangent, last_tangent);
  leftBound_ = referenceLine_.toFrenetPolyline(left_boundary_pts);
  rightBound_ = referenceLine_.toFrenetPolyline(right_boundary_pts);
}

BoundaryDistances Corridor::signedDistancesAt(
    const RealType arc_length) const noexcept {
  return std::make_pair(leftBound_.deviationAt(arc_length),
                        rightBound_.deviationAt(arc_length));
}

RealType Corridor::widthAt(const RealType arc_length) const noexcept {
  return leftBound_.deviationAt(arc_length) +
         std::abs(rightBound_.deviationAt(arc_length));
}

RealType Corridor::centerOffset(const RealType arc_length) const noexcept {
  const BoundaryDistances distances = signedDistancesAt(arc_length);
  return (distances.first + distances.second) * 0.5;
}

RealType Corridor::lengthReferenceLine() const noexcept {
  return referenceLine_.GetTotalLength();
}

FrenetFrame2D Corridor::FrenetFrame(
    const CartesianPoint2D& position) const noexcept {
  return referenceLine_.getFrenetPositionWithFrame(position).frame;
}

FrenetPositionWithFrame Corridor::getFrenetPositionWithFrame(
    const CartesianPoint2D& position) const noexcept {
  return referenceLine_.getFrenetPositionWithFrame(position);
}

std::ostream& operator<<(std::ostream& os, const CorridorPtr& corridor) {
  using namespace std;
  os << "Corridor " << corridor->id() << "\n";
  os << corridor->referenceLine_ << "\n";
  os << corridor->leftBound_ << "\n";
  os << corridor->rightBound_ << "\n";
  // os << "Attributes:";
  // for (const auto& o : corridor.attributes()) {
  //   os << " [" << o.first << ": " << o.second << "]";
  // }
  return os;
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