#include <boost/python.hpp>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"
#include "corridor/frenet_types.h"

// Python API wrapper
#include "corridor_assignment_wrapper.hpp"
#include "corridor_wrapper.hpp"
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
  // Corridor Wrapper
  // ///////////////////////////////////////////////////////////////////////////
  py::class_<CorridorWrapper>("CorridorWrapper",
                              py::init<int, const py::list&, const py::list&>())
      .def("get_polylines", &CorridorWrapper::GetCartesianPolylinesLines)
      .def("to_frenet_state_vector", &CorridorWrapper::ToFrenetStateVector)
      .def("to_frenet_state", &CorridorWrapper::ToFrenetState);

  py::def("TestCorridorHandle", &TestCorridorHandle);

  // ///////////////////////////////////////////////////////////////////////////
  //  UT Wrapper
  // ///////////////////////////////////////////////////////////////////////////

  py::class_<FlatCartesianPositionAndCovMat2D>(
      "FlatCartesianPositionAndCovMat2D")
      .def_readwrite("x", &FlatCartesianPositionAndCovMat2D::x)
      .def_readwrite("y", &FlatCartesianPositionAndCovMat2D::y)
      .def_readwrite("var_x", &FlatCartesianPositionAndCovMat2D::var_x)
      .def_readwrite("var_y", &FlatCartesianPositionAndCovMat2D::var_y)
      .def_readwrite("cov_xy", &FlatCartesianPositionAndCovMat2D::cov_xy);

  py::class_<FlatPolarPositionAndCovMat2D>("FlatPolarPositionAndCovMat2D")
      .def_readwrite("r", &FlatPolarPositionAndCovMat2D::r)
      .def_readwrite("phi", &FlatPolarPositionAndCovMat2D::phi)
      .def_readwrite("var_r", &FlatPolarPositionAndCovMat2D::var_r)
      .def_readwrite("var_phi", &FlatPolarPositionAndCovMat2D::var_phi)
      .def_readwrite("cov_rphi", &FlatPolarPositionAndCovMat2D::cov_rphi);

  def("cartesian_to_polar_2d", &pyCartesianToPolarTransformation2D);
  def("polar_to_cartesian_2d", &pyPolarToCartesianTransformation2D);
  boost::python::def("ut_cartesian_to_polar_2d",
                     &UnscentedTransformationPolarCoordinate2D);

  py::class_<FlatCartesianStateAndCovMat2D>("FlatCartesianStateAndCovMat2D")
      .def_readwrite("x", &FlatCartesianStateAndCovMat2D::x)
      .def_readwrite("y", &FlatCartesianStateAndCovMat2D::y)
      .def_readwrite("vx", &FlatCartesianStateAndCovMat2D::vx)
      .def_readwrite("vy", &FlatCartesianStateAndCovMat2D::vy)
      .def_readwrite("var_x", &FlatCartesianStateAndCovMat2D::var_x)
      .def_readwrite("var_y", &FlatCartesianStateAndCovMat2D::var_y)
      .def_readwrite("var_vx", &FlatCartesianStateAndCovMat2D::var_vx)
      .def_readwrite("var_vy", &FlatCartesianStateAndCovMat2D::var_vy)
      .def_readwrite("cov_xy", &FlatCartesianStateAndCovMat2D::cov_xy)
      .def_readwrite("cov_xvx", &FlatCartesianStateAndCovMat2D::cov_xvx)
      .def_readwrite("cov_xvy", &FlatCartesianStateAndCovMat2D::cov_xvy)
      .def_readwrite("cov_yvx", &FlatCartesianStateAndCovMat2D::cov_yvx)
      .def_readwrite("cov_yvy", &FlatCartesianStateAndCovMat2D::cov_yvy)
      .def_readwrite("cov_vxvy", &FlatCartesianStateAndCovMat2D::cov_vxvy);

  py::class_<FlatFrenetStateAndCovMat2D>("FlatFrenetStateAndCovMat2D")
      .def_readwrite("l", &FlatFrenetStateAndCovMat2D::l)
      .def_readwrite("d", &FlatFrenetStateAndCovMat2D::d)
      .def_readwrite("vl", &FlatFrenetStateAndCovMat2D::vl)
      .def_readwrite("vd", &FlatFrenetStateAndCovMat2D::vd)
      .def_readwrite("var_d", &FlatFrenetStateAndCovMat2D::var_d)
      .def_readwrite("var_l", &FlatFrenetStateAndCovMat2D::var_l)
      .def_readwrite("var_vd", &FlatFrenetStateAndCovMat2D::var_vd)
      .def_readwrite("var_vl", &FlatFrenetStateAndCovMat2D::var_vl)
      .def_readwrite("cov_ld", &FlatFrenetStateAndCovMat2D::cov_ld)
      .def_readwrite("cov_lvl", &FlatFrenetStateAndCovMat2D::cov_lvl)
      .def_readwrite("cov_lvd", &FlatFrenetStateAndCovMat2D::cov_lvd)
      .def_readwrite("cov_dvl", &FlatFrenetStateAndCovMat2D::cov_dvl)
      .def_readwrite("cov_dvd", &FlatFrenetStateAndCovMat2D::cov_dvd)
      .def_readwrite("cov_vlvd", &FlatFrenetStateAndCovMat2D::cov_vlvd);

  // UnscentedStateTransformation
  boost::python::def("ut_cartesian_frenet_transformation",
                     &UnscentedStateTransformation);

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
