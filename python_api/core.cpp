#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_interpolation_2d.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/cubic_spline/cubic_spline_utilities.h"
#include "corridor/frenet_types.h"
#include "unscented_transformation_wrapper.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
namespace cr = corridor;
namespace cs = cr::cubic_spline;
namespace ut = cr::unscented_transformation;

// TODO: replace this by later exposing the spline class to python
cs::DataMatrix<cr::RealType> spline_data_matrix_;

template <typename T>
inline std::vector<T> to_std_vector(const py::object& iterable) {
  return std::vector<T>(py::stl_input_iterator<T>(iterable),
                        py::stl_input_iterator<T>());
}

inline py::dict to_py_dict(const cs::SplineCoefficients2d& coeffs) {
  py::dict dict;

  py::list a_x, b_x, c_x, d_x;
  py::list a_y, b_y, c_y, d_y;
  for (const auto& c : coeffs) {
    a_x.append(c.a_x);
    b_x.append(c.b_x);
    c_x.append(c.c_x);
    d_x.append(c.d_x);

    a_y.append(c.a_y);
    b_y.append(c.b_y);
    c_y.append(c.c_y);
    d_y.append(c.d_y);
  }
  dict["a_x"] = a_x;
  dict["b_x"] = b_x;
  dict["c_x"] = c_x;
  dict["d_x"] = d_x;

  dict["a_y"] = a_y;
  dict["b_y"] = b_y;
  dict["c_y"] = c_y;
  dict["d_y"] = d_y;

  return dict;
}

inline py::dict to_py_dict(const cr::CartesianPoint2D& vector2d) {
  py::dict py_dict;
  py_dict["x"] = vector2d.x();
  py_dict["y"] = vector2d.y();
  return py_dict;
}

inline py::list to_py_list(const cr::FrenetFrames2D& frenet_frames) {
  py::list py_frenet_frames;
  for (const auto& ff : frenet_frames) {
    py::dict py_frenet_frame;
    py_frenet_frame["segm_index"] = ff.frenet_base().segment_info.idx;
    py_frenet_frame["segm_arc_length"] =
        ff.frenet_base().segment_info.relative_arc_length;
    py_frenet_frame["origin"] = to_py_dict(ff.origin());
    py_frenet_frame["tangent"] = to_py_dict(ff.tangent());
    py_frenet_frame["normal"] = to_py_dict(ff.normal());
    py_frenet_frames.append(py_frenet_frame);
  }

  return py_frenet_frames;
}

py::dict CreateSplineParams(const py::list& node_x, const py::list& node_y) {
  using namespace corridor;
  std::vector<RealType> x_vec = to_std_vector<RealType>(node_x);
  std::vector<RealType> y_vec = to_std_vector<RealType>(node_y);

  // Create natural spline from points
  cr::CartesianPoints2D points;
  for (int i = 0, n = x_vec.size(); i < n; i++) {
    points.emplace_back(x_vec[i], y_vec[i]);
  }
  spline_data_matrix_ = cs::naturalSplineDataMatrixFromPoints(points);

  py::list arc_lengths;
  for (cs::DataIdx i = 0, max_idx = spline_data_matrix_.cols(); i < max_idx;
       i++) {
    arc_lengths.append(
        spline_data_matrix_.col(i)[cs::BasicDataTypes::kArcLength]);
  }

  // Construct spline coefficients
  cs::SplineCoefficients2d coefficients =
      cs::SplineCoefficientsFromDataMatrix(spline_data_matrix_);

  py::dict py_dict = to_py_dict(coefficients);
  py_dict["arc_length"] = arc_lengths;

  return py_dict;
}

py::list ConstructFrenetFrames(const cr::RealType target_x,
                               const cr::RealType target_y) {
  using namespace corridor;

  if (spline_data_matrix_.cols() == 0) {
    std::cout << "Data matrix not yet set! Call 'CreateSplineParams' first."
              << std::endl;
    return py::list();
  }

  CartesianPoint2D target_point;
  target_point << target_x, target_y;

  const FrenetFrames2D frenet_frames =
      cs::ConstructFrenetFrames(spline_data_matrix_, target_point);

  return to_py_list(frenet_frames);
}

// /////////////////////////////////////////////////////////////////////////////
// Uncertainty Transformation
// /////////////////////////////////////////////////////////////////////////////

py::list CartesianToPolarTransformation2D(const cr::RealType x1,
                                          const cr::RealType x2) {
  Eigen::Vector2d polar_vector = ut::CartesianToPolarTransformation2D({x1, x2});
  py::list polar_coordinates;
  polar_coordinates.append(polar_vector(0));
  polar_coordinates.append(polar_vector(1));
  return polar_coordinates;
}

py::list PolarToCartesianTransformation2D(const cr::RealType radius,
                                          const cr::RealType phi) {
  Eigen::Vector2d cartesian_vector =
      ut::PolarToCartesianTransformation2D({radius, phi});
  py::list cartesian_coordinates;
  cartesian_coordinates.append(cartesian_vector(0));
  cartesian_coordinates.append(cartesian_vector(1));
  return cartesian_coordinates;
}

ut::FlatPolarStateAndCovMat2D UnscentedTransformationPolarCoordinate2D(
    const ut::FlatCartesianStateAndCovMat2D& cartesian_state) {
  cr::StateMeanAndCovarianceMatrix polar_state =
      ut::UnscentedTransformationPolarCoordinates(ut::Convert(cartesian_state));
  return ut::Convert(polar_state);
}

BOOST_PYTHON_MODULE(PYTHON_API_MODULE_NAME) {  // NOLINT
  def("create_spline_params", &CreateSplineParams);
  def("construct_frenet_frames", &ConstructFrenetFrames);

  def("cartesian_to_polar_2d", &CartesianToPolarTransformation2D);
  def("polar_to_cartesian_2d", &PolarToCartesianTransformation2D);
  boost::python::def("ut_cartesian_to_polar_2d",
                     &UnscentedTransformationPolarCoordinate2D);

  py::class_<ut::FlatCartesianStateAndCovMat2D>("FlatCartesianStateAndCovMat2D")
      .def_readwrite("x", &ut::FlatCartesianStateAndCovMat2D::x)
      .def_readwrite("y", &ut::FlatCartesianStateAndCovMat2D::y)
      .def_readwrite("var_x", &ut::FlatCartesianStateAndCovMat2D::var_x)
      .def_readwrite("var_y", &ut::FlatCartesianStateAndCovMat2D::var_y)
      .def_readwrite("cov_xy", &ut::FlatCartesianStateAndCovMat2D::cov_xy);

  py::class_<ut::FlatPolarStateAndCovMat2D>("FlatPolarStateAndCovMat2D")
      .def_readwrite("r", &ut::FlatPolarStateAndCovMat2D::r)
      .def_readwrite("phi", &ut::FlatPolarStateAndCovMat2D::phi)
      .def_readwrite("var_r", &ut::FlatPolarStateAndCovMat2D::var_r)
      .def_readwrite("var_phi", &ut::FlatPolarStateAndCovMat2D::var_phi)
      .def_readwrite("cov_rphi", &ut::FlatPolarStateAndCovMat2D::cov_rphi);
}
