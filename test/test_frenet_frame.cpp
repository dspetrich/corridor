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
  EXPECT_EQ(frenet_frame.frenet_base().id, InvalidId);
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

TEST_F(FrenetFrameTest, testFrenetStateConversion) {  // NOLINT

  const CartesianPoint2D position(5.0, 4.0);
  const CartesianVector2D velocity(3.0, 5.0);
  const CovarianceMatrix2D cm_position(1.0, 1.5);
  const CovarianceMatrix2D cm_velocity(2.0, 3.0);

  CartesianState2D cartesian_state(position, velocity, cm_position,
                                   cm_velocity);
  FrenetState2D frenet_state =
      frenet_frame_.FromCartesianState(cartesian_state);

  EXPECT_FLOAT_EQ(frenet_state.position().l(), 15.34);
  EXPECT_FLOAT_EQ(frenet_state.position().d(), 1.0);
  EXPECT_FLOAT_EQ(frenet_state.velocity().l(), 3.0);
  EXPECT_FLOAT_EQ(frenet_state.velocity().d(), 5.0);

  // Since the Frenet frame shares the same axes as the cartesian frame, the
  // covariances are identical
  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().position().ll(),
                  cm_position.xx());
  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().position().dd(),
                  cm_position.yy());
  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().position().ld(),
                  cm_position.xy());

  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().velocity().ll(),
                  cm_velocity.xx());
  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().velocity().dd(),
                  cm_velocity.yy());
  EXPECT_FLOAT_EQ(frenet_state.covarianceMatrix().velocity().ld(),
                  cm_velocity.xy());

  const auto abs_velocity = frenet_state.abs_velocity();
  const auto orientation = frenet_state.orientation();

  EXPECT_FLOAT_EQ(abs_velocity.value, 6.02515);
  EXPECT_FLOAT_EQ(orientation.value, 1.0169548);

  EXPECT_FLOAT_EQ(abs_velocity.variance(), 2.8107188);
  EXPECT_FLOAT_EQ(orientation.variance(), 0.0669813);
}