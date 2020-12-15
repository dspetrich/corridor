#include <gtest/gtest.h>

#include "corridor/cubic_spline/cubic_interpolation_2d.h"

using namespace corridor;
using namespace cubic_spline;

class Interpolation2dTest : public ::testing::Test {
 protected:
  void SetUp() override {
    refline_.emplace_back(1, 0);
    refline_.emplace_back(2, 3);
    refline_.emplace_back(3, 6);
    refline_.emplace_back(4, 5);
  }

 public:
  CartesianPoints2D refline_;
};

TEST_F(Interpolation2dTest, Interpolation) {
  DataMatrix<RealType> data = NaturalSplineDataMatrixFromPoints(refline_);

  // Test end point of evaluation with nodes
  for (DataIdx i = 0, max_idx = data.cols() - 1; i < max_idx; i++) {
    auto l = data.col(i + 1)[kArcLength] - data.col(i)[kArcLength];
    auto p = Coefficients2d(data.col(i), data.col(i + 1)).evaluatePosition(l);
    auto expect_p = data.block<2, 1>(kPoint_x, i + 1);
    EXPECT_LE((p - expect_p).norm(), 1e-6);
  }
}
