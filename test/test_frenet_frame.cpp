#include <gtest/gtest.h>

#include "corridor/frenet_types.h"

using namespace corridor;
using namespace cubic_spline;

class FrenetFrameTest : public ::testing::Test {
  void SetUp() override {
    frenet_base_ = FrenetBase2D(23, 12.34, 0.0, -.01, 0.002, {3, 2.34});
    origin_ = CartesianPoint2D(2.0, 3.0);
    tangent_ = CartesianVector2D(1.0, 0.0);
    normal_ = CartesianVector2D(0.0, 1.0);

    frenet_frame_ = FrenetFrame2D(frenet_base_, origin_, tangent_, normal_);
  }

 public:
  FrenetBase2D frenet_base_;
  CartesianPoint2D origin_;
  CartesianVector2D tangent_;
  CartesianVector2D normal_;

  FrenetFrame2D frenet_frame_;
};

TEST_F(FrenetFrameTest, EmptyConstructor) {  // NOLINT
  FrenetFrame2D frenet_frame;
  EXPECT_EQ(frenet_frame.frenet_base().id, InvalId);
}

TEST_F(FrenetFrameTest, Constructor) {  // NOLINT
  EXPECT_EQ(frenet_frame_.frenet_base().id, frenet_base_.id);
  EXPECT_EQ(frenet_frame_.frenet_base().segment_info.idx,
            frenet_base_.segment_info.idx);
  EXPECT_FLOAT_EQ(frenet_frame_.frenet_base().arc_length,
                  frenet_base_.arc_length);
  EXPECT_FLOAT_EQ(frenet_frame_.frenet_base().segment_info.relative_arc_length,
                  frenet_base_.segment_info.relative_arc_length);

  CartesianPoint2D delta_origin = frenet_frame_.origin() - origin_;
  EXPECT_FLOAT_EQ(delta_origin.norm(), 0.f);

  CartesianVector2D delta_tangent = frenet_frame_.tangent() - tangent_;
  EXPECT_FLOAT_EQ(delta_tangent.norm(), 0.f);

  CartesianVector2D delta_normal = frenet_frame_.normal() - normal_;
  EXPECT_FLOAT_EQ(delta_normal.norm(), 0.f);

  EXPECT_FLOAT_EQ(frenet_frame_.curvature(), -0.01);
}

TEST_F(FrenetFrameTest, Conversion) {  // NOLINT
  // Check point/position conversion
  CartesianPoint2D cartesian_point(3.f, 2.f);
  FrenetPoint2D frenet_point =
      frenet_frame_.FromCartesianPoint(cartesian_point);

  CartesianPoint2D cartesian_point2 =
      frenet_frame_.FromFrenetPoint(frenet_point);

  CartesianPoint2D delta = cartesian_point - cartesian_point2;

  EXPECT_FLOAT_EQ(delta.norm(), 0.f);
  EXPECT_FLOAT_EQ(frenet_point.l(), 13.34);
  EXPECT_FLOAT_EQ(frenet_point.d(), -1.0);

  // Check vector conversion
  CartesianVector2D cartesian_vector(3.f, 2.f);
  FrenetVector2D frenet_vector =
      frenet_frame_.FromCartesianVector(cartesian_vector);
  CartesianVector2D cartesian_vector2 =
      frenet_frame_.FromFrenetVector(frenet_vector);

  CartesianVector2D delta2 = cartesian_vector - cartesian_vector2;
  EXPECT_FLOAT_EQ(delta.norm(), 0.f);

  EXPECT_FLOAT_EQ(frenet_vector.l(), 3.0);
  EXPECT_FLOAT_EQ(frenet_vector.d(), 2.0);
}