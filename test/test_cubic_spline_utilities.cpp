#include <gtest/gtest.h>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"

using namespace corridor;
using namespace cubic_spline;

class CubicSplineUtilitiesTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // heart shaped spline
    refline_.emplace_back(0, -6);
    refline_.emplace_back(5, 0);
    refline_.emplace_back(4, 3);
    refline_.emplace_back(2, 4);
    refline_.emplace_back(1, 3);
    refline_.emplace_back(0, 0);
    refline_.emplace_back(-1, 3);
    refline_.emplace_back(-2, 4);
    refline_.emplace_back(-4, 3);
    refline_.emplace_back(-5, 0);
    refline_.emplace_back(0, -6);

    point_ << 0.f, -2.f;
  }

 public:
  CartesianPoints2D refline_;
  CartesianPoint2D point_;
};

TEST_F(CubicSplineUtilitiesTest, FindSegmentCandidates) {  // NOLINT
  DataMatrix<RealType> data = naturalSplineDataMatrixFromPoints(refline_);

  SegmentInfoVector<DataIdx, RealType> segment_points;
  const bool matched = FindSegmentCandidates(data, point_, &segment_points);
  ASSERT_TRUE(matched);
  EXPECT_EQ(segment_points.size(), 3);
  EXPECT_FLOAT_EQ(segment_points[0].idx, 0);
  EXPECT_FLOAT_EQ(segment_points[1].idx, 4);
  EXPECT_FLOAT_EQ(segment_points[2].idx, 9);
}

TEST_F(CubicSplineUtilitiesTest, FrenetFrame) {  // NOLINT
  DataMatrix<RealType> data = naturalSplineDataMatrixFromPoints(refline_);

  const FrenetFrames2D frenet_frames = ConstructFrenetFrames(data, point_);

  EXPECT_EQ(frenet_frames.size(), 3);
  const auto& ff1 = frenet_frames[0];
  const auto& ff2 = frenet_frames[1];
  const auto& ff3 = frenet_frames[2];

  EXPECT_EQ(ff1.frenet_base().segment_info.idx, 0);
  EXPECT_EQ(ff2.frenet_base().segment_info.idx, 4);
  EXPECT_EQ(ff3.frenet_base().segment_info.idx, 9);

  EXPECT_FLOAT_EQ(ff1.frenet_base().arc_length, 2.3020294);
  EXPECT_FLOAT_EQ(ff2.frenet_base().arc_length, 18.290928);
  EXPECT_FLOAT_EQ(ff3.frenet_base().arc_length, 34.279823);

  // Origin
  CartesianPoint2D expected;
  expected << 2.01695, -4.47608;
  EXPECT_NEAR((expected - ff1.origin()).norm(), 0.f, 1e-5);
  expected << 0.f, 0.f;
  EXPECT_NEAR((expected - ff2.origin()).norm(), 0.f, 1e-5);
  expected << -2.01695, -4.47608;
  EXPECT_NEAR((expected - ff3.origin()).norm(), 0.f, 1e-5);

  // tangent
  expected << 0.775401, 0.631469;
  EXPECT_NEAR((expected - ff1.tangent()).norm(), 0.f, 1e-5);
  expected << -1.f, 0.f;
  EXPECT_NEAR((expected - ff2.tangent()).norm(), 0.f, 1e-5);
  expected << 0.775401, -0.631469;
  EXPECT_NEAR((expected - ff3.tangent()).norm(), 0.f, 1e-5);

  // normal
  expected << -0.631469, 0.775401;
  EXPECT_NEAR((expected - ff1.normal()).norm(), 0.f, 1e-5);
  expected << 0.f, -1.f;
  EXPECT_NEAR((expected - ff2.normal()).norm(), 0.f, 1e-5);
  expected << 0.631469, 0.775401;
  EXPECT_NEAR((expected - ff3.normal()).norm(), 0.f, 1e-5);
}
