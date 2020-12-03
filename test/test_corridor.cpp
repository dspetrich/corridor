#include <gtest/gtest.h>

#include <cmath>

#include "corridor/corridor.h"
#include "lanelet_mock.hpp"

using namespace corridor;

class CorridorTest : public ::testing::Test {
 protected:
  void SetUp() override {
    id_ = 123;
    expected_arclength_ = M_SQRT2 * 6;
    straight_corridor_ = Corridor(id_, straight_lanelet_.centerline,
                                  straight_lanelet_.left_boundary,
                                  straight_lanelet_.right_boundary);
  }

 public:
  StraightLanelet straight_lanelet_;
  FirstStraightLaneletSection first_lanelet_;
  SecondStraightLaneletSection second_lanelet_;

  Corridor straight_corridor_;

  IdType id_;
  RealType expected_arclength_;
};

TEST_F(CorridorTest, createEmptyCorridor) {
  Corridor corridor;
  EXPECT_EQ(corridor.id(), InvalidId);
  EXPECT_EQ(corridor.lengthReferenceLine(), 0.0);
}

TEST_F(CorridorTest, createStraightCorridor) {
  EXPECT_EQ(straight_corridor_.id(), id_);
  EXPECT_FLOAT_EQ(straight_corridor_.lengthReferenceLine(),
                  expected_arclength_);

  EXPECT_FLOAT_EQ(straight_corridor_.widthAt(0), 2 * M_SQRT2);
  EXPECT_FLOAT_EQ(straight_corridor_.widthAt(expected_arclength_ / 2.0),
                  2 * M_SQRT2);
  EXPECT_FLOAT_EQ(straight_corridor_.widthAt(expected_arclength_), 2 * M_SQRT2);

  EXPECT_FLOAT_EQ(
      straight_corridor_.signedDistancesAt(expected_arclength_ / 3.0).first,
      M_SQRT2);
  EXPECT_FLOAT_EQ(
      straight_corridor_.signedDistancesAt(expected_arclength_ / 3.0).second,
      -M_SQRT2);

  EXPECT_FLOAT_EQ(straight_corridor_.centerOffset(0), 0.0);
  EXPECT_FLOAT_EQ(straight_corridor_.centerOffset(expected_arclength_ / 2.0),
                  0.0);
  EXPECT_FLOAT_EQ(straight_corridor_.centerOffset(expected_arclength_), 0.0);

  EXPECT_FLOAT_EQ(straight_corridor_.curvatureAt(0), 0.0);
  EXPECT_FLOAT_EQ(straight_corridor_.curvatureAt(expected_arclength_ / 2.0),
                  0.0);
  EXPECT_FLOAT_EQ(straight_corridor_.curvatureAt(expected_arclength_), 0.0);
}

TEST_F(CorridorTest, frenetStateTransformation) {
  CartesianPoint2D position(4.0, 3.0);
  FrenetFrame2D frenet_frame = straight_corridor_.FrenetFrame(position);

  EXPECT_EQ(frenet_frame.frenet_base().id, id_);
  EXPECT_NEAR(frenet_frame.frenet_base().arc_length, 3.5 * M_SQRT2, 1e-3);
  EXPECT_FLOAT_EQ(frenet_frame.frenet_base().orientation, M_PI / 4.0);
  EXPECT_FLOAT_EQ(frenet_frame.frenet_base().curvature, 0.0);
  EXPECT_FLOAT_EQ(frenet_frame.frenet_base().curvature_change_rate, 0.0);
  EXPECT_EQ(frenet_frame.frenet_base().segment_info.idx, 1);
  EXPECT_NEAR(frenet_frame.frenet_base().segment_info.relative_arc_length,
              0.5 * M_SQRT2, 1e-3);
}

TEST_F(CorridorTest, testCorridorSequence) {
  const auto first_corridor_ptr = std::make_shared<Corridor>(
      1, first_lanelet_.centerline, first_lanelet_.left_boundary,
      first_lanelet_.right_boundary);
  const auto second_corridor_ptr = std::make_shared<Corridor>(
      2, second_lanelet_.centerline, second_lanelet_.left_boundary,
      second_lanelet_.right_boundary);

  CorridorSequence corridor_sequence({first_corridor_ptr, second_corridor_ptr});

  EXPECT_FALSE(corridor_sequence.empty());
  EXPECT_EQ(corridor_sequence.size(), 2);
  EXPECT_EQ(corridor_sequence.totalLength(), 20.0);

  EXPECT_FLOAT_EQ(corridor_sequence.signedDistancesAt(18).first, 2.0);
  EXPECT_FLOAT_EQ(corridor_sequence.signedDistancesAt(18).second, -2.0);
  EXPECT_FLOAT_EQ(corridor_sequence.widthAt(18), 4.0);

  // Test FrenetFrame conversion with different start arc_length which results
  // in different start corridors within the corridor sequence
  const CartesianPoint2D position(2, 0.4);

  FrenetPositionWithFrame frenet_data1 =
      corridor_sequence.getFrenetPositionWithFrame(position, -100.0);
  EXPECT_EQ(frenet_data1.frame.frenet_base().id, 1);
  EXPECT_EQ(frenet_data1.frame.frenet_base().segment_info.idx, 0);

  FrenetPositionWithFrame frenet_data2 =
      corridor_sequence.getFrenetPositionWithFrame(position);
  EXPECT_EQ(frenet_data2.frame.frenet_base().id, 1);
  EXPECT_EQ(frenet_data2.frame.frenet_base().segment_info.idx, 0);

  FrenetPositionWithFrame frenet_data3 =
      corridor_sequence.getFrenetPositionWithFrame(position, 40);
  EXPECT_EQ(frenet_data3.frame.frenet_base().id, 1);
  EXPECT_EQ(frenet_data3.frame.frenet_base().segment_info.idx, 0);

  // Last corridor in the sequence
  const CartesianPoint2D position2(16.0, -1.7);
  FrenetPositionWithFrame frenet_data5 =
      corridor_sequence.getFrenetPositionWithFrame(position2, 100.0);
  EXPECT_EQ(frenet_data5.frame.frenet_base().id, 2);
  EXPECT_EQ(frenet_data5.frame.frenet_base().segment_info.idx, 1);

  FrenetPositionWithFrame frenet_data6 =
      corridor_sequence.getFrenetPositionWithFrame(position2, 3.0);
  EXPECT_EQ(frenet_data6.frame.frenet_base().id, 2);
  EXPECT_EQ(frenet_data6.frame.frenet_base().segment_info.idx, 1);

  FrenetPositionWithFrame frenet_data7 =
      corridor_sequence.getFrenetPositionWithFrame(position2, -100.0);
  EXPECT_EQ(frenet_data7.frame.frenet_base().id, 2);
  EXPECT_EQ(frenet_data7.frame.frenet_base().segment_info.idx, 1);
}