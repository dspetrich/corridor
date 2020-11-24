#include <gtest/gtest.h>

#include "corridor/cubic_spline/cubic_spline.h"
#include "corridor/frenet_types.h"

using namespace corridor;
using namespace cubic_spline;

class CubicSplineTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // heart shaped spline
    x_vec_ = {0, 5, 4, 2, 1, 0, -1.f, -2.f, -4.f, -5.f, 0.f};
    y_vec_ = {-6, 0, 3, 4, 3, 0, 3.f, 4.f, 3.f, 0.f, -6.f};
  }

 public:
  std::vector<RealType> x_vec_;
  std::vector<RealType> y_vec_;
};

TEST_F(CubicSplineTest, constructor) {  // NOLINT
  CubicSpline default_constructor_spline;
  EXPECT_EQ(InvalId, default_constructor_spline.GetId());
  EXPECT_EQ(0, default_constructor_spline.GetSize());
  EXPECT_FLOAT_EQ(0.f, default_constructor_spline.GetTotalLength());

  // set id
  CubicSpline cubic_spline_2(100);
  EXPECT_EQ(100, cubic_spline_2.GetId());
  EXPECT_EQ(0, cubic_spline_2.GetSize());
  EXPECT_FLOAT_EQ(0.f, cubic_spline_2.GetTotalLength());

  // set id and epsilon
  std::size_t epsilon = 20;
  CubicSpline cubic_spline_3(200, epsilon);
  EXPECT_EQ(200, cubic_spline_3.GetId());
  EXPECT_EQ(0, cubic_spline_3.GetSize());
  EXPECT_FLOAT_EQ(0.f, cubic_spline_3.GetTotalLength());

  // invalid (empty) points
  std::vector<RealType> x_vec, y_vec;
  CubicSpline cubic_spline_4(300, x_vec, y_vec);  // invlad points
  EXPECT_EQ(InvalId, cubic_spline_4.GetId());
  EXPECT_EQ(0, cubic_spline_4.GetSize());
  EXPECT_FLOAT_EQ(0.f, cubic_spline_4.GetTotalLength());
}

TEST_F(CubicSplineTest, Evaluate) {  // NOLINT
  IdType id = 23;
  CubicSpline cubic_spline(id, x_vec_, y_vec_);
  EXPECT_EQ(id, cubic_spline.GetId());
  EXPECT_EQ(11, cubic_spline.GetSize());
  EXPECT_FLOAT_EQ(36.581856, cubic_spline.GetTotalLength());
}

TEST_F(CubicSplineTest, FrenetFrame) {  // NOLINT
  IdType id = 23;
  CubicSpline cubic_spline(id, x_vec_, y_vec_);

  CartesianPoint2D point;
  point << 0.f, -2.f;
  FrenetFrames2D frenet_frames = cubic_spline.FrenetFrames(point);

  EXPECT_EQ(frenet_frames.size(), 3);
}
