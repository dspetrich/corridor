#include <boost/python.hpp>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"
#include "corridor/frenet_types.h"

// Python API wrapper
#include "corridor_assignment_wrapper.hpp"
#include "cubic_spline_wrapper.hpp"
#include "unscented_transformation_wrapper.hpp"

namespace py = boost::python;

BOOST_PYTHON_MODULE(PYTHON_API_MODULE_NAME) {  // NOLINT

  // ///////////////////////////////////////////////////////////////////////////
  // Cubic Spline
  // ///////////////////////////////////////////////////////////////////////////
  py::class_<CubicSpline>("CubicSpline")
      .def("naturalSplineParameter", &CubicSpline::naturalSplineParameter)
      .def("clampedSplineParameter", &CubicSpline::clampedSplineParameter)
      .def("constructFrenetFrames", &CubicSpline::constructFrenetFrames);

  // ///////////////////////////////////////////////////////////////////////////
  // Polar Uncertainty Transformation Wrapper
  // ///////////////////////////////////////////////////////////////////////////

  py::class_<FlatCartesianStateAndCovMat2D>("FlatCartesianStateAndCovMat2D")
      .def_readwrite("x", &FlatCartesianStateAndCovMat2D::x)
      .def_readwrite("y", &FlatCartesianStateAndCovMat2D::y)
      .def_readwrite("var_x", &FlatCartesianStateAndCovMat2D::var_x)
      .def_readwrite("var_y", &FlatCartesianStateAndCovMat2D::var_y)
      .def_readwrite("cov_xy", &FlatCartesianStateAndCovMat2D::cov_xy);

  py::class_<FlatPolarStateAndCovMat2D>("FlatPolarStateAndCovMat2D")
      .def_readwrite("r", &FlatPolarStateAndCovMat2D::r)
      .def_readwrite("phi", &FlatPolarStateAndCovMat2D::phi)
      .def_readwrite("var_r", &FlatPolarStateAndCovMat2D::var_r)
      .def_readwrite("var_phi", &FlatPolarStateAndCovMat2D::var_phi)
      .def_readwrite("cov_rphi", &FlatPolarStateAndCovMat2D::cov_rphi);

  def("cartesian_to_polar_2d", &pyCartesianToPolarTransformation2D);
  def("polar_to_cartesian_2d", &pyPolarToCartesianTransformation2D);
  boost::python::def("ut_cartesian_to_polar_2d",
                     &UnscentedTransformationPolarCoordinate2D);

  // ///////////////////////////////////////////////////////////////////////////
  // Corridor Assignment Wrapper
  // ///////////////////////////////////////////////////////////////////////////
  py::class_<FlatCorridorRelatedFeatures>("CorridorAssignmentFeature")
      .def_readwrite("l", &FlatCorridorRelatedFeatures::l)
      .def_readwrite("d", &FlatCorridorRelatedFeatures::d)
      .def_readwrite("sigma_l", &FlatCorridorRelatedFeatures::sigma_l)
      .def_readwrite("sigma_d", &FlatCorridorRelatedFeatures::sigma_d)
      .def_readwrite("obj_width_ratio",
                     &FlatCorridorRelatedFeatures::obj_width_ratio)
      .def_readwrite("obj_length_ratio",
                     &FlatCorridorRelatedFeatures::obj_length_ratio)
      .def_readwrite("corridor_width",
                     &FlatCorridorRelatedFeatures::corridor_width)
      .def_readwrite("corridor_length",
                     &FlatCorridorRelatedFeatures::corridor_length);

  //! callable functions in Python
  py::def("LateralAssignmentConfidence", &LateralAssignmentConfidence);
  py::def("LongitudinalAssignmentConfidence",
          &LongitudinalAssignmentConfidence);

  py::def("evaluateIntegralLineWidthGaussian",
          &evaluateIntegralLineWidthGaussian);

  py::def("MovingConfidence", &MovingConfidence);
  py::def("RelativeDirectionConfidence", &RelativeDirectionConfidence);
}
