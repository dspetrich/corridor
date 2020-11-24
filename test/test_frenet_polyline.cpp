#include <gtest/gtest.h>

#include <cmath>

#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/frenet_types.h"

using namespace corridor;
using namespace cubic_spline;

struct LaneletMock {
  CartesianPoints2D centerline;
  CartesianPoints2D left_boundary;
  CartesianPoints2D right_boundary;
};

struct StraightLanelet : public LaneletMock {
  void fillCenterline() {
    centerline.emplace_back(0.f, 0.f);
    centerline.emplace_back(3.f, 3.f);
    centerline.emplace_back(6.f, 6.f);
  }
  void fillLeftBoundary() {
    left_boundary.emplace_back(-2.f, 0.f);
    left_boundary.emplace_back(-1.f, 1.f);
    left_boundary.emplace_back(0.f, 2.f);
    left_boundary.emplace_back(1.f, 3.f);
    left_boundary.emplace_back(2.f, 4.f);
    left_boundary.emplace_back(3.f, 5.f);
    left_boundary.emplace_back(4.f, 6.f);
    left_boundary.emplace_back(5.f, 7.f);
    left_boundary.emplace_back(6.f, 8.f);
  }
  void fillRightBoundary() {
    right_boundary.emplace_back(0.f, -2.f);
    right_boundary.emplace_back(1.f, -1.f);
    right_boundary.emplace_back(2.f, 0.f);
    right_boundary.emplace_back(3.f, 1.f);
    right_boundary.emplace_back(4.f, 2.f);
    right_boundary.emplace_back(5.f, 3.f);
    right_boundary.emplace_back(6.f, 4.f);
    right_boundary.emplace_back(7.f, 5.f);
    right_boundary.emplace_back(8.f, 6.f);
  }
  StraightLanelet() {
    fillCenterline();
    fillLeftBoundary();
    fillRightBoundary();
  }
};

struct CurvedLanelet : public LaneletMock {
  void fillCenterline() {
    centerline.emplace_back(0.f, 0.f);
    centerline.emplace_back(4.f, 5.f);
    centerline.emplace_back(8.f, 0.f);
  }
  void fillLeftBoundary() {
    left_boundary.emplace_back(-2.f, -1.f);
    left_boundary.emplace_back(-1.f, 1.f);
    left_boundary.emplace_back(0.f, 3.f);
    left_boundary.emplace_back(1.f, 5.f);
    left_boundary.emplace_back(4.f, 6.f);
    left_boundary.emplace_back(7.f, 5.f);
    left_boundary.emplace_back(8.f, 3.f);
    left_boundary.emplace_back(9.f, 1.f);
    left_boundary.emplace_back(10.f, -1.f);
  }
  void fillRightBoundary() {
    right_boundary.emplace_back(0.f, -3.f);
    right_boundary.emplace_back(1.f, -1.f);
    right_boundary.emplace_back(2.f, 1.f);
    right_boundary.emplace_back(3.f, 3.f);
    right_boundary.emplace_back(4.f, 4.f);
    right_boundary.emplace_back(5.f, 3.f);
    right_boundary.emplace_back(6.f, 1.f);
    right_boundary.emplace_back(7.f, -1.f);
    right_boundary.emplace_back(8.f, -3.f);
  }
  CurvedLanelet() {
    fillCenterline();
    fillLeftBoundary();
    fillRightBoundary();
  }
};

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
  EXPECT_FLOAT_EQ(left_polyline[1].l(), 0.31910014);
  EXPECT_FLOAT_EQ(left_polyline[2].l(), 2.0599406);
  EXPECT_FLOAT_EQ(left_polyline[3].l(), 3.8812499);
  EXPECT_FLOAT_EQ(left_polyline[4].l(), 6.605114);
  EXPECT_FLOAT_EQ(left_polyline[5].l(), 9.3289785);
  EXPECT_FLOAT_EQ(left_polyline[6].l(), 11.150288);
  EXPECT_FLOAT_EQ(left_polyline[7].l(), 12.891129);
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
  EXPECT_FLOAT_EQ(right_polyline[2].l(), 1.4562994);
  EXPECT_FLOAT_EQ(right_polyline[3].l(), 3.5452135);
  EXPECT_FLOAT_EQ(right_polyline[4].l(), 6.6175189);
  EXPECT_FLOAT_EQ(right_polyline[5].l(), 9.6650152);
  EXPECT_FLOAT_EQ(right_polyline[6].l(), 11.753929);
  EXPECT_FLOAT_EQ(right_polyline[7].l(), 13.621993);
  EXPECT_FLOAT_EQ(right_polyline[8].l(), 15.857287);

  EXPECT_FLOAT_EQ(right_polyline[0].d(), -1.41176474);
  EXPECT_FLOAT_EQ(right_polyline[1].d(), -1.35294116);
  EXPECT_FLOAT_EQ(right_polyline[2].d(), -1.28179097);
  EXPECT_FLOAT_EQ(right_polyline[3].d(), -1.06581819);
  EXPECT_FLOAT_EQ(right_polyline[4].d(), -1.0000018);
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
  EXPECT_FLOAT_EQ(d_10, 1.6338902);

  auto right_polyline =
      curved_cubic_spline_.toFrenetPolyline(curved_lanelet_.right_boundary);
  std::cout << right_polyline.deviationAt(10.0) << std::endl;

  const auto d2_tooShort = right_polyline.deviationAt(-3.0);
  EXPECT_FLOAT_EQ(d2_tooShort, right_polyline[0].d());

  const auto d2_tooLarge = right_polyline.deviationAt(16.0);
  EXPECT_FLOAT_EQ(d2_tooLarge, right_polyline[8].d());

  const auto d2_10 = right_polyline.deviationAt(10.0);
  EXPECT_LT(d2_10, right_polyline[5].d());
  EXPECT_GT(d2_10, right_polyline[6].d());
  EXPECT_FLOAT_EQ(d2_10, -1.1004523);
}