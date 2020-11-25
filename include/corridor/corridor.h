#pragma once

#include <lanelet2_core/primitives/Lanelet.h>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/frenet_types.h"

namespace corridor {

enum class CorridorType { kNone = 0, kLaneSegment = 1, kCrosswalk = 2 };

class Corridor;
using CorridorPtr = std::shared_ptr<Corridor>;
using CorridorPtrs = std::vector<CorridorPtr>;

/**
 * @brief Distance of the left and right boundary of a corridor to its
 * reference line.
 * first: distance to left bound
 * second: distance to right bound
 */
using BoundaryDistances = std::pair<RealType, RealType>;

class Corridor {
 public:
  Corridor() : id_(InvalidId), type_(CorridorType::kNone) {}
  // Corridor(const ll::ConstLanelet& lanelet);

  //! Get the unique id of underlying lane let
  IdType id() const noexcept { return id_; }

  BoundaryDistances signedDistancesAt(const RealType arc_length) const noexcept;
  RealType widthAt(const RealType arc_length) const noexcept;
  RealType centerOffset(const RealType arc_length) const noexcept;
  RealType curvatureAt(const RealType arc_length) const noexcept {
    return referenceLine_.GetCurvatureAt(arc_length);
  }
  RealType lengthReferenceLine() const noexcept;

  FrenetFrame2D FrenetFrame(const CartesianPoint2D& position) const noexcept;

  FrenetPositionWithFrame getFrenetPositionWithFrame(
      const CartesianPoint2D& position) const noexcept;

  // Introspection
  friend std::ostream& operator<<(std::ostream& os,
                                  const CorridorPtr& corridor);

 private:
  // Id and type of the corridor
  IdType id_;
  CorridorType type_;

  // Reference line for the frenet frame. Not necessarily a centerline, but has
  // to be located between the left and right boundary.
  cubic_spline::CubicSpline referenceLine_;

  // Left and right side of the corridor are independent from the sampling of
  // the reference line and independent from each other to assure maximal
  // flexibility
  FrenetPolyline leftBound_;
  FrenetPolyline rightBound_;
};

// Introspection
std::ostream& operator<<(std::ostream& os, const CorridorPtr& corridor);

/**
 * @brief an ordered sequence of corridors
 *
 * A CorridorSequence is a collection of Corridors that behaves like
 * like one single Corridor.
 */
class CorridorSequence {
 public:
  CorridorSequence(const CorridorPtrs& corridor_ptrs = CorridorPtrs()) {
    RealType offset = 0.0;
    std::transform(corridor_ptrs.begin(), corridor_ptrs.end(),
                   std::inserter(corridor_map_, corridor_map_.end()),
                   [&offset](const CorridorPtr c) {
                     RealType prev_offset = offset;
                     offset += c->lengthReferenceLine();
                     return std::make_pair(prev_offset, c);
                   });
  }

  //! Returns wether this holds any corridors
  bool empty() const noexcept { return corridor_map_.empty(); }

  //! Returns number of corridors
  size_t size() const noexcept { return corridor_map_.size(); }

  BoundaryDistances signedDistancesAt(const RealType arc_length) const noexcept;
  RealType widthAt(const RealType arc_length) const noexcept;
  RealType centerOffsetAt(const RealType arc_length) const noexcept;
  RealType curvatureAt(const RealType arc_length) const noexcept;
  RealType totalLength() const noexcept;

  FrenetPositionWithFrame getFrenetPositionWithFrame(
      const CartesianPoint2D& position,
      const RealType start_arc_length = 0.0) const noexcept;

  // Introspection
  friend std::ostream& operator<<(std::ostream& os,
                                  const CorridorPtr& corridor);

 private:
  // key: offset of the corridor to the first corridor
  // value: pointer to the corridor object
  using CorridorSequenceMap = std::map<RealType, CorridorPtr>;
  CorridorSequenceMap corridor_map_;

  // Get corridor with contains the arc_length.
  // If it is smaller than zero it will return the first corridor, if it is
  // bigger than the totalLength() it will return the last CorridorPtr
  CorridorSequenceMap::const_iterator get(const RealType arc_length) const {
    CorridorSequenceMap::const_iterator it =
        corridor_map_.lower_bound(arc_length);
    assert(!corridor_map_.empty());
    if (it == corridor_map_.end()) {
      // If key is larger than the biggest key, return last map entry
      return std::prev(it);
    }
    return it;
  }

  FrenetPositionWithFrame getFrenetPositionWithFrame(
      const CartesianPoint2D& position,
      CorridorSequenceMap::const_iterator iter) const;
};

using CorridorPath = CorridorPtrs;
using CorridorPaths = std::vector<CorridorPath>;

std::ostream& operator<<(std::ostream& os, const CorridorPath& path);
std::ostream& operator<<(std::ostream& os, const CorridorPaths& paths);

}  // namespace corridor
