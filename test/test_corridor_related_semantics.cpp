#include <gtest/gtest.h>

#include <cmath>

#include "corridor/corridor_assignment/corridor_related_semantics.h"

using namespace corridor;

TEST(CorridorRelatedSemantic, vectorInitialization) {
  SemanticLabelSet semantic_labels({SemanticLabel::kDownstream});

  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kDownstream));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kUndefined));

  EXPECT_FALSE(semantic_labels.hasLabel(SemanticLabel::kUpstream));
  EXPECT_FALSE(semantic_labels.hasLabel(SemanticLabel::kTowardsLeft));
  EXPECT_FALSE(semantic_labels.hasLabel(SemanticLabel::kTowardsRight));

  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kUndefined), 1.0);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kDownstream), 0.0);
}

TEST(CorridorRelatedSemantic, SemanticLabelPairs) {
  SemanticLabelPairs slp{{SemanticLabel::kUpstream, 0.2},
                         {SemanticLabel::kTowardsLeft, 0.1},
                         {SemanticLabel::kTowardsRight, 0.1},
                         {SemanticLabel::kDownstream, 0.1}};

  SemanticLabelSet semantic_labels(slp);

  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kUndefined));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kUpstream));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kDownstream));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kTowardsLeft));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kTowardsRight));

  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kUndefined), 0.5);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kUpstream), 0.2);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kDownstream), 0.1);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kTowardsLeft), 0.1);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kTowardsRight), 0.1);
}

TEST(CorridorRelatedSemantic, SemanticLabelPairsAboveOne) {
  SemanticLabelPairs slp{{SemanticLabel::kUpstream, 0.8},
                         {SemanticLabel::kTowardsLeft, 0.4},
                         {SemanticLabel::kTowardsRight, 0.2},
                         {SemanticLabel::kDownstream, 0.1}};

  SemanticLabelSet semantic_labels(slp);

  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kUndefined));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kUpstream));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kDownstream));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kTowardsLeft));
  EXPECT_TRUE(semantic_labels.hasLabel(SemanticLabel::kTowardsRight));

  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kUndefined), 0.0);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kUpstream), 0.53333336);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kDownstream), 0.06666667);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kTowardsLeft), 0.26666668);
  EXPECT_FLOAT_EQ(semantic_labels.at(SemanticLabel::kTowardsRight), 0.13333334);
}