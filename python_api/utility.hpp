#pragma once

#include <boost/python.hpp>

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline_coefficients.h"
#include "corridor/frenet_types.h"

namespace py = boost::python;
namespace cs = corridor::cubic_spline;

// /////////////////////////////////////////////////////////////////////////////
// Utility functions
// /////////////////////////////////////////////////////////////////////////////

template <typename T>
inline std::vector<T> to_std_vector(const py::object& iterable) {
  return std::vector<T>(py::stl_input_iterator<T>(iterable),
                        py::stl_input_iterator<T>());
}

inline corridor::CartesianVector2D to_cartesian_vector(
    const py::list& py_list) {
  using namespace corridor;
  std::vector<RealType> std_vec = to_std_vector<RealType>(py_list);
  return CartesianVector2D(std_vec.front(), std_vec.back());
}

inline corridor::CartesianPoints2D toCartesianPoints(const py::list& node_x,
                                                     const py::list& node_y) {
  using namespace corridor;
  std::vector<RealType> x_vec = to_std_vector<RealType>(node_x);
  std::vector<RealType> y_vec = to_std_vector<RealType>(node_y);

  corridor::CartesianPoints2D points;
  for (int i = 0, n = x_vec.size(); i < n; i++) {
    points.emplace_back(x_vec[i], y_vec[i]);
  }
  return points;
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

inline py::dict to_py_dict(const corridor::CartesianPoint2D& vector2d) {
  py::dict py_dict;
  py_dict["x"] = vector2d.x();
  py_dict["y"] = vector2d.y();
  return py_dict;
}

inline py::list to_py_list(const corridor::FrenetFrames2D& frenet_frames) {
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

inline std::pair<py::list, py::list> to_py_lists(
    const corridor::CartesianPoints2D& points) {
  py::list point_x;
  py::list point_y;
  for (const corridor::CartesianPoint2D& p : points) {
    point_x.append(p.x());
    point_y.append(p.y());
  }
  return std::make_pair(point_x, point_y);
}

corridor::CartesianStateVector2D convert(
    const py::list& cartesian_state_vector) {
  using namespace corridor;
  return CartesianStateVector2D(
      py::extract<RealType>(cartesian_state_vector[0]),
      py::extract<RealType>(cartesian_state_vector[1]),
      py::extract<RealType>(cartesian_state_vector[2]),
      py::extract<RealType>(cartesian_state_vector[3]));
}

py::list convert(const corridor::FrenetStateVector2D& frenet_state_vector) {
  py::list py_frenet_state_vector;
  py_frenet_state_vector.append(frenet_state_vector.l());
  py_frenet_state_vector.append(frenet_state_vector.d());
  py_frenet_state_vector.append(frenet_state_vector.vl());
  py_frenet_state_vector.append(frenet_state_vector.vd());
  return py_frenet_state_vector;
}

struct FlatCartesianPositionAndCovMat2D {
  corridor::RealType x;
  corridor::RealType y;
  corridor::RealType var_x;
  corridor::RealType var_y;
  corridor::RealType cov_xy;
};

struct FlatCartesianStateAndCovMat2D {
  corridor::RealType x;
  corridor::RealType y;
  corridor::RealType vx;
  corridor::RealType vy;

  corridor::RealType var_x;
  corridor::RealType var_y;
  corridor::RealType var_vx;
  corridor::RealType var_vy;

  corridor::RealType cov_xy;
  corridor::RealType cov_xvx;
  corridor::RealType cov_xvy;
  corridor::RealType cov_yvx;
  corridor::RealType cov_yvy;
  corridor::RealType cov_vxvy;
};

struct FlatPolarPositionAndCovMat2D {
  corridor::RealType r;
  corridor::RealType phi;
  corridor::RealType var_r;
  corridor::RealType var_phi;
  corridor::RealType cov_rphi;
};

struct FlatFrenetStateAndCovMat2D {
  corridor::RealType l;
  corridor::RealType d;
  corridor::RealType vl;
  corridor::RealType vd;

  corridor::RealType var_l;
  corridor::RealType var_d;
  corridor::RealType var_vl;
  corridor::RealType var_vd;

  corridor::RealType cov_ld;
  corridor::RealType cov_lvl;
  corridor::RealType cov_lvd;
  corridor::RealType cov_dvl;
  corridor::RealType cov_dvd;
  corridor::RealType cov_vlvd;
};

corridor::StateMeanAndCovarianceMatrix Convert(
    const FlatCartesianPositionAndCovMat2D& flat_cartesian_state) {
  corridor::StateMeanAndCovarianceMatrix state(2);
  state.mean << flat_cartesian_state.x, flat_cartesian_state.y;
  state.covMat(0, 0) = flat_cartesian_state.var_x;
  state.covMat(1, 1) = flat_cartesian_state.var_y;
  state.covMat(1, 0) = flat_cartesian_state.cov_xy;
  state.covMat(0, 1) = flat_cartesian_state.cov_xy;
  return state;
}

corridor::CartesianState2D Convert(
    const FlatCartesianStateAndCovMat2D& flat_cartesian_state) {
  using namespace corridor;
  const CartesianPoint2D position(flat_cartesian_state.x,
                                  flat_cartesian_state.y);
  const CartesianVector2D velocity(flat_cartesian_state.vx,
                                   flat_cartesian_state.vy);
  const CovarianceMatrix2D cm_position(flat_cartesian_state.var_x,
                                       flat_cartesian_state.var_y,
                                       flat_cartesian_state.cov_xy);
  const CovarianceMatrix2D cm_velocity(flat_cartesian_state.var_vx,
                                       flat_cartesian_state.var_vy,
                                       flat_cartesian_state.cov_vxvy);
  const CovarianceMatrix2D cm_pos_vel(
      flat_cartesian_state.cov_xvx, flat_cartesian_state.cov_yvy,
      flat_cartesian_state.cov_xvy, flat_cartesian_state.cov_yvx);
  return CartesianState2D(position, velocity, cm_position, cm_velocity,
                          cm_pos_vel);
}

FlatFrenetStateAndCovMat2D Convert(
    const corridor::FrenetState2D& frenet_state) {
  FlatFrenetStateAndCovMat2D flat_frenet_state;
  flat_frenet_state.l = frenet_state.l();
  flat_frenet_state.d = frenet_state.d();
  flat_frenet_state.vl = frenet_state.vl();
  flat_frenet_state.vd = frenet_state.vd();

  const corridor::FrenetStateCovarianceMatrix2D& cov_mat =
      frenet_state.covarianceMatrix();

  flat_frenet_state.var_l = cov_mat.ll();
  flat_frenet_state.var_d = cov_mat.dd();
  flat_frenet_state.cov_ld = cov_mat.ld();

  flat_frenet_state.var_vl = cov_mat.vlvl();
  flat_frenet_state.var_vd = cov_mat.vdvd();
  flat_frenet_state.cov_vlvd = cov_mat.vlvd();

  flat_frenet_state.cov_lvl = cov_mat.lvl();
  flat_frenet_state.cov_lvd = cov_mat.lvd();
  flat_frenet_state.cov_dvl = cov_mat.dvl();
  flat_frenet_state.cov_dvd = cov_mat.dvd();

  return flat_frenet_state;
}

FlatPolarPositionAndCovMat2D Convert(
    const corridor::StateMeanAndCovarianceMatrix& polar_state) {
  FlatPolarPositionAndCovMat2D flat_polar_state;
  flat_polar_state.r = polar_state.mean(0);
  flat_polar_state.phi = polar_state.mean(1);
  flat_polar_state.var_r = polar_state.covMat(0, 0);
  flat_polar_state.var_phi = polar_state.covMat(1, 1);
  flat_polar_state.cov_rphi = polar_state.covMat(1, 0);
  return flat_polar_state;
};
