#pragma once

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "corridor/basic_types.h"

namespace corridor {

// /////////////////////////////////////////////////////////////////////////////
// Unordered map utilities
// /////////////////////////////////////////////////////////////////////////////

struct EnumClassHash {
  template <typename T>
  std::size_t operator()(T t) const {
    return static_cast<std::size_t>(t);
  }
};

struct AddValues {
  template <class Value, class Pair>
  Value operator()(Value value, const Pair &pair) const {
    return value + pair.second;
  }
};

// /////////////////////////////////////////////////////////////////////////////
// Semantic labels
// /////////////////////////////////////////////////////////////////////////////

// More labels possible
enum class SemanticLabel : uint64_t {
  kUndefined = 0,  //!< No valid label
  kDownstream,     //!< Towards moving direction of corridor
  kUpstream,       //!< Against moving direction of corridor
  kTowardsLeft,    //!< Towards left side
  kTowardsRight,   //!< Towards right side
};
using SemanticLabels = std::vector<SemanticLabel>;

using SemanticLabelPair = std::pair<SemanticLabel, RealType>;
using SemanticLabelPairs = std::vector<SemanticLabelPair>;

/**
 * @brief Container which associates a Semantic Label with its confidence value
 *        key: SemanticLabel
 *        value: confidence in [0, 1]
 */
using SemanticLabelMap =
    std::unordered_map<SemanticLabel, RealType, EnumClassHash>;

/**
 * @brief Container which holds a specific set of semantic labels and their
 * confidence value.
 * This container represents discrete probability distribution over all
 * included labels. Meaning, that the sum over all confidences is guaranteed to
 * be one.
 */
class SemanticLabelSet {
 public:
  SemanticLabelSet(const SemanticLabels &labels) {
    // Initialize label map with provided keys
    std::transform(labels.begin(), labels.end(),
                   std::inserter(label_map_, label_map_.end()),
                   [](SemanticLabel l) { return std::make_pair(l, 0.0); });
    // Set initial value to kUndefined, no matter if it was provided in the
    // given labels or not
    label_map_[SemanticLabel::kUndefined] = 1.0;
  }
  SemanticLabelSet(const SemanticLabelPairs &label_pairs) {
    // Initialize label map with provided key-value pairs
    label_map_ = SemanticLabelMap(label_pairs.begin(), label_pairs.end());
    initializeNormalization();
  }

  /**
   * @brief Set confidence value of a specific label if present in the label
   * map. Return false when the label is not present and doesn't change the map.
   *
   * @param[in] label: label which will be altered
   * @param[in] confidence: new confidence value for the label
   * @return true: label is present in the map and was altered
   * @return false: label is not present in the map, nothing was changed
   */
  bool setLabel(const SemanticLabel label, const RealType confidence) {
    auto label_iter = label_map_.find(label);
    if (label_iter == label_map_.end()) {
      return false;
    }
    label_iter->second = confidence;
    return true;
  }

  /**
   * @brief Check if a certain label is present in this map
   *
   * @param label
   * @return true: label is present
   * @return false: label is not present
   */
  bool hasLabel(const SemanticLabel label) const {
    const auto label_iter = label_map_.find(label);
    if (label_iter != label_map_.end()) {
      return true;
    }
    return false;
  }

  RealType at(const SemanticLabel label) const { return label_map_.at(label); }

  void normalize() {
    RealType sum =
        std::accumulate(label_map_.begin(), label_map_.end(), 0.0, AddValues());
    for (auto &p : label_map_) {
      p.second /= sum;
    }
  }

 private:
  SemanticLabelMap label_map_;

  void initializeNormalization() {
    auto undefined_label_iter = label_map_.find(SemanticLabel::kUndefined);

    if (undefined_label_iter != label_map_.end()) {
      // Undefined is present, so just normalization is missing
      normalize();
      return;
    }

    // Add this point we need to add the kUndefined label and normalize.
    const RealType sum =
        std::accumulate(label_map_.begin(), label_map_.end(), 0.0, AddValues());

    if (sum >= 1.0) {
      // Sum of all labels are already sufficient, so normalization is
      // sufficient
      label_map_[SemanticLabel::kUndefined] = 0.0;
      for (auto &p : label_map_) {
        p.second /= sum;
      }
      return;
    }

    // If the sum is below 1, add the remaining confidence to kUndefined.
    label_map_[SemanticLabel::kUndefined] = 1 - sum;
  }
};

}  // namespace corridor