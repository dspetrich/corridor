#include <gtest/gtest.h>

#include "corridor/unscented_transformation/sigma_points.h"
#include "corridor/unscented_transformation/unscented_transformation.h"

using namespace corridor;
using namespace unscented_transformation;

class UnscentedTransformationTest : public ::testing::Test {
 protected:
  void SetUp() override {
    initial_state_ = StateMeanAndCovarianceMatrix(3);
    sigma_points_ = MerweScaledSigmaPoints(3);
    initial_state_.mean << 5, 2.2, M_PI_2;
    initial_state_.covMat(0, 0) = 0.5;
    initial_state_.covMat(1, 1) = 0.1;
    initial_state_.covMat(2, 2) = M_PI / 8.0;
  }

 public:
  StateMeanAndCovarianceMatrix initial_state_;
  MerweScaledSigmaPoints sigma_points_;
};

TEST_F(UnscentedTransformationTest, initialization) {
  EXPECT_EQ(initial_state_.mean.size(), 3);
  EXPECT_FLOAT_EQ(initial_state_.mean.norm(), 5.68396);
  EXPECT_EQ(initial_state_.covMat.size(), 9);
  EXPECT_FLOAT_EQ(initial_state_.covMat.determinant(), 0.019634955);
}

TEST_F(UnscentedTransformationTest, SigmaPoints) {
  const auto& sigmas = sigma_points_.generateSigmaPoints(initial_state_.mean,
                                                         initial_state_.covMat);
  EXPECT_EQ((sigmas.col(0) - initial_state_.mean).norm(), 0.0);
  EXPECT_EQ(sigmas.cols(), 7);
  EXPECT_FLOAT_EQ(sigma_points_.weightsMean().sum(), 1.0);
  // Merwe weights don't necessarily sum to one!
  EXPECT_FLOAT_EQ(sigma_points_.weightsCovMat().sum(), 3.99);
}

TEST_F(UnscentedTransformationTest, UT_boxTransformation) {
  const auto& sigmas = sigma_points_.generateSigmaPoints(initial_state_.mean,
                                                         initial_state_.covMat);

  Eigen::MatrixXd transformed_sigmas(initial_state_.mean.size(), sigmas.cols());
  for (int i = 0; i < sigmas.cols(); i++) {
    BoxDimension box(sigmas(0, i), sigmas(1, i));
    const auto projection = box.projection(sigmas(2, i));
    transformed_sigmas(0, i) = projection.first;
    transformed_sigmas(1, i) = projection.second;
    transformed_sigmas(2, i) = sigmas(2, i);
  }

  StateMeanAndCovarianceMatrix transformed_state =
      EstimateStateMeanAndCovarianceMatrix(transformed_sigmas,
                                           sigma_points_.weightsMean(),
                                           sigma_points_.weightsCovMat(), 2);

  std::cout << "/* transformed_state */" << std::endl;
  std::cout << transformed_state.mean << std::endl;

  std::cout << "/* transformed_covMat */" << std::endl;
  std::cout << transformed_state.covMat << std::endl;

  // Nice transformation function, like the projection of a box to the frame
  // axes
  //
}