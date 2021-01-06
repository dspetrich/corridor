#include <gtest/gtest.h>

#include <cmath>

#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/frenet_types.h"
#include "lanelet_mock.hpp"

using namespace corridor;
using namespace cubic_spline;

class FrenetPolylineTest : public ::testing::Test {
  void SetUp() override {
    straight_cubic_spline_ = CubicSpline(1, straight_lanelet_.centerline);
    curved_cubic_spline_ = CubicSpline(2, curved_lanelet_.centerline);
  }

 public:
  StraightLanelet straight_lanelet_;
  CurvedLanelet curved_lanelet_;
  CubicSpline straight_cubic_spline_;
  CubicSpline curved_cubic_spline_;
};

TEST_F(FrenetPolylineTest, Straight) {
  auto left_polyline =
      straight_cubic_spline_.toFrenetPolyline(straight_lanelet_.left_boundary);

  EXPECT_EQ(left_polyline.size(), straight_lanelet_.left_boundary.size());
  for (int i = 0, n = left_polyline.size(); i < n; i++) {
    EXPECT_NEAR(left_polyline[i].l(), M_SQRT2 * (i - 1), 1e-12);
    EXPECT_NEAR(left_polyline[i].d(), M_SQRT2, 1e-12);
  }

  auto right_polyline =
      straight_cubic_spline_.toFrenetPolyline(straight_lanelet_.right_boundary);
  EXPECT_EQ(right_polyline.size(), straight_lanelet_.left_boundary.size());
  for (int i = 0, n = right_polyline.size(); i < n; i++) {
    EXPECT_NEAR(right_polyline[i].l(), M_SQRT2 * (i - 1), 1e-12);
    EXPECT_NEAR(right_polyline[i].d(), -M_SQRT2, 1e-12);
  }
}

TEST_F(FrenetPolylineTest, Curved) {
  auto left_polyline =
      curved_cubic_spline_.toFrenetPolyline(curved_lanelet_.left_boundary);

  EXPECT_EQ(left_polyline.size(), curved_lanelet_.left_boundary.size());

  EXPECT_FLOAT_EQ(left_polyline[0].l(), -1.82352948);
  EXPECT_FLOAT_EQ(left_polyline[1].l(), 0.31897849);
  EXPECT_FLOAT_EQ(left_polyline[2].l(), 2.0597572);
  EXPECT_FLOAT_EQ(left_polyline[3].l(), 3.8812857);
  EXPECT_FLOAT_EQ(left_polyline[4].l(), 6.605114);
  EXPECT_FLOAT_EQ(left_polyline[5].l(), 9.3288383);
  EXPECT_FLOAT_EQ(left_polyline[6].l(), 11.150427);
  EXPECT_FLOAT_EQ(left_polyline[7].l(), 12.891085);
  EXPECT_FLOAT_EQ(left_polyline[8].l(), 15.033757);

  EXPECT_FLOAT_EQ(left_polyline[0].d(), 1.29411769);
  EXPECT_FLOAT_EQ(left_polyline[1].d(), 1.35307455);
  EXPECT_FLOAT_EQ(left_polyline[2].d(), 1.44882464);
  EXPECT_FLOAT_EQ(left_polyline[3].d(), 1.74184847);
  EXPECT_FLOAT_EQ(left_polyline[4].d(), 1.f);
  EXPECT_FLOAT_EQ(left_polyline[5].d(), 1.74184883);
  EXPECT_FLOAT_EQ(left_polyline[6].d(), 1.44882488);
  EXPECT_FLOAT_EQ(left_polyline[7].d(), 1.35307467);
  EXPECT_FLOAT_EQ(left_polyline[8].d(), 1.29411793);

  auto right_polyline =
      curved_cubic_spline_.toFrenetPolyline(curved_lanelet_.right_boundary);

  EXPECT_EQ(right_polyline.size(), curved_lanelet_.right_boundary.size());

  EXPECT_FLOAT_EQ(right_polyline[0].l(), -2.64705896);
  EXPECT_FLOAT_EQ(right_polyline[1].l(), -0.411764711);
  EXPECT_FLOAT_EQ(right_polyline[2].l(), 1.4561728);
  EXPECT_FLOAT_EQ(right_polyline[3].l(), 3.5453079);
  EXPECT_FLOAT_EQ(right_polyline[4].l(), 6.6058898);
  EXPECT_FLOAT_EQ(right_polyline[5].l(), 9.6649895);
  EXPECT_FLOAT_EQ(right_polyline[6].l(), 11.753894);
  EXPECT_FLOAT_EQ(right_polyline[7].l(), 13.621993);
  EXPECT_FLOAT_EQ(right_polyline[8].l(), 15.857287);

  EXPECT_FLOAT_EQ(right_polyline[0].d(), -1.41176474);
  EXPECT_FLOAT_EQ(right_polyline[1].d(), -1.35294116);
  EXPECT_FLOAT_EQ(right_polyline[2].d(), -1.28179097);
  EXPECT_FLOAT_EQ(right_polyline[3].d(), -1.06581819);
  EXPECT_FLOAT_EQ(right_polyline[4].d(), -1.0);
  EXPECT_FLOAT_EQ(right_polyline[5].d(), -1.06581807);
  EXPECT_FLOAT_EQ(right_polyline[6].d(), -1.28179049);
  EXPECT_FLOAT_EQ(right_polyline[7].d(), -1.35294092);
  EXPECT_FLOAT_EQ(right_polyline[8].d(), -1.4117645);
}

TEST_F(FrenetPolylineTest, Interpolation) {
  auto left_polyline =
      curved_cubic_spline_.toFrenetPolyline(curved_lanelet_.left_boundary);

  const auto d_tooShort = left_polyline.deviationAt(-2.0);
  EXPECT_FLOAT_EQ(d_tooShort, left_polyline[0].d());

  const auto d_tooLarge = left_polyline.deviationAt(16.0);
  EXPECT_FLOAT_EQ(d_tooLarge, left_polyline[8].d());

  const auto d_10 = left_polyline.deviationAt(10.0);
  EXPECT_LT(d_10, left_polyline[5].d());
  EXPECT_GT(d_10, left_polyline[6].d());
  EXPECT_FLOAT_EQ(d_10, 1.6338842);

  auto right_polyline =
      curved_cubic_spline_.toFrenetPolyline(curved_lanelet_.right_boundary);

  const auto d2_tooShort = right_polyline.deviationAt(-3.0);
  EXPECT_FLOAT_EQ(d2_tooShort, right_polyline[0].d());

  const auto d2_tooLarge = right_polyline.deviationAt(16.0);
  EXPECT_FLOAT_EQ(d2_tooLarge, right_polyline[8].d());

  const auto d2_10 = right_polyline.deviationAt(10.0);
  EXPECT_LT(d2_10, right_polyline[5].d());
  EXPECT_GT(d2_10, right_polyline[6].d());
  EXPECT_FLOAT_EQ(d2_10, -1.1004552);
}