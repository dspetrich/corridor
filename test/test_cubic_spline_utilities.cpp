#include <gtest/gtest.h>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_segment_root_finding.h"
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
  DataMatrix<RealType> data = NaturalSplineDataMatrixFromPoints(refline_);

  SegmentInfoVector<DataIdx, RealType> segment_points;
  const bool matched = FindSegmentCandidates(data, point_, &segment_points);
  ASSERT_TRUE(matched);
  EXPECT_EQ(segment_points.size(), 3);
  EXPECT_FLOAT_EQ(segment_points[0].idx, 0);
  EXPECT_FLOAT_EQ(segment_points[1].idx, 4);
  EXPECT_FLOAT_EQ(segment_points[2].idx, 9);
}

TEST_F(CubicSplineUtilitiesTest, FindProjectionArcLength) {  // NOLINT
  DataMatrix<RealType> data = NaturalSplineDataMatrixFromPoints(refline_);

  SegmentInfoVector<DataIdx, RealType> segment_candidates;
  const bool matched = FindSegmentCandidates(data, point_, &segment_candidates);
  ASSERT_TRUE(matched);
  EXPECT_EQ(segment_candidates.size(), 3);
  EXPECT_FLOAT_EQ(segment_candidates[0].idx, 0);
  EXPECT_FLOAT_EQ(segment_candidates[1].idx, 4);
  EXPECT_FLOAT_EQ(segment_candidates[2].idx, 9);

  // Calculate data segments
  const DataSegment<RealType>& data_segment_0 =
      data.block<kSize, 2>(kPoint_x, segment_candidates[0].idx);
  const Coefficients2d segment_coeffs_0(data_segment_0.col(0),
                                        data_segment_0.col(1));

  // Find pependicular projection of point onto segments
  for (auto& segment : segment_candidates) {
    const DataSegment<RealType>& data_segment =
        data.block<kSize, 2>(kPoint_x, segment.idx);
    const Coefficients2d segment_coeffs(data_segment.col(0),
                                        data_segment.col(1));
    if (matched) {
      const auto bisection =
          BisectionMethod(segment_coeffs, segment.relative_arc_length, point_);
      const auto brent =
          BrentsMethod(segment_coeffs, segment.relative_arc_length, point_);

      if (segment.idx == 4) {
        // Here, the point is perpendicular to the endpoint of the segment.
        // Therefore Brent's Methode doesn't find the root, since the sign of
        // function values doesn't change. This is expected and is handled by
        // LimitArcLengthToSegmentLimits, which prevents the call of Brent's
        // Method in the regular processing.
        ASSERT_FALSE(brent.first);
        continue;
      }

      ASSERT_TRUE(bisection.first);
      ASSERT_TRUE(brent.first);

      EXPECT_NEAR(bisection.second, brent.second, 2e-3);
    }
  }
}

TEST_F(CubicSplineUtilitiesTest, FrenetFrame) {  // NOLINT
  DataMatrix<RealType> data = NaturalSplineDataMatrixFromPoints(refline_);

  const FrenetFrames2D frenet_frames = ConstructFrenetFrames(data, point_);

  EXPECT_EQ(frenet_frames.size(), 3);
  const auto& ff1 = frenet_frames[0];
  const auto& ff2 = frenet_frames[1];
  const auto& ff3 = frenet_frames[2];

  EXPECT_EQ(ff1.frenet_base().segment_info.idx, 0);
  EXPECT_EQ(ff2.frenet_base().segment_info.idx, 4);
  EXPECT_EQ(ff3.frenet_base().segment_info.idx, 9);

  EXPECT_FLOAT_EQ(ff1.frenet_base().arc_length, 2.302582);
  EXPECT_FLOAT_EQ(ff2.frenet_base().arc_length, 18.290928);
  EXPECT_FLOAT_EQ(ff3.frenet_base().arc_length, 34.280846);

  // Origin
  CartesianPoint2D expected;
  expected << 2.01741, -4.4757;
  EXPECT_NEAR((expected - ff1.origin()).norm(), 0.f, 1e-5);
  expected << 0.f, 0.f;
  EXPECT_NEAR((expected - ff2.origin()).norm(), 0.f, 1e-5);
  expected << -2.0161, -4.47676;
  EXPECT_NEAR((expected - ff3.origin()).norm(), 0.f, 1e-5);

  // tangent
  expected << 0.775383, 0.631491;
  EXPECT_NEAR((expected - ff1.tangent()).norm(), 0.f, 1e-5);
  expected << -1.f, 0.f;
  EXPECT_NEAR((expected - ff2.tangent()).norm(), 0.f, 1e-5);
  expected << 0.775431, -0.631432;
  EXPECT_NEAR((expected - ff3.tangent()).norm(), 0.f, 1e-5);

  // normal
  expected << -0.631491, 0.775383;
  EXPECT_NEAR((expected - ff1.normal()).norm(), 0.f, 1e-5);
  expected << 0.f, -1.f;
  EXPECT_NEAR((expected - ff2.normal()).norm(), 0.f, 1e-5);
  expected << 0.631432, 0.775431;
  EXPECT_NEAR((expected - ff3.normal()).norm(), 0.f, 1e-5);
}
