#pragma once

#include "corridor/basic_types.h"

namespace corridor {
namespace unscented_transformation {

Eigen::Vector2d CartesianToPolarTransformation2D(
    const Eigen::Vector2d& cartesian_vector);

Eigen::Vector2d PolarToCartesianTransformation2D(
    const Eigen::Vector2d& polar_vector);

void ToPolarCoordinates2D(const Eigen::Vector2d& initial_x,
                          const Eigen::Matrix2d& initial_P,
                          PolarVector2D* resulting_x,
                          PolarCovarianceMatrix2D* resulting_P);

void FromPolarCoordinates2D(const PolarVector2D& initial_x,
                            const PolarCovarianceMatrix2D& initial_P,
                            Eigen::Vector2d* resulting_x,
                            Eigen::Vector2d* resulting_P);

}  // namespace unscented_transformation
}  // namespace corridor