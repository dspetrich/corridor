#pragma once

#include <iostream>  // std::cout

#include "corridor/basic_types.h"
#include "corridor/cartesian_types.h"
#include "corridor/cubic_spline/cubic_spline_types.h"

namespace corridor {

struct FrenetPoint2D : public BasicPoint2D {
  FrenetPoint2D(void) : BasicPoint2D() {}
  FrenetPoint2D(const RealType l, const RealType d) : BasicPoint2D(l, d) {}

  // This constructor allows you to construct FrenetPoint2D from Eigen
  // expressions
  template <typename OtherDerived>
  FrenetPoint2D(const Eigen::MatrixBase<OtherDerived>& other)
      : BasicPoint2D(other) {}

  typedef Eigen::Matrix<RealType, 2, 1, Eigen::DontAlign> Base;

  // This method allows you to assign Eigen expressions to FrenetPoint2D
  template <typename OtherDerived>
  FrenetPoint2D& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Easy access to the frenet coordinates
  const RealType l() const { return (*this)[0]; }
  const RealType d() const { return (*this)[1]; }

  RealType& l() { return (*this)[0]; }
  RealType& d() { return (*this)[1]; }

  const RealType d_value() const { return std::abs((*this)[0]); }
};
using FrenetVector2D = FrenetPoint2D;
using FrenetPoints2D = std::vector<FrenetPoint2D>;
using FrenetVectors2D = std::vector<FrenetVector2D>;

// /////////////////////////////////////////////////////////////////////////////
// Frenet Base, Frenet Frame and Frenet Polyline
// /////////////////////////////////////////////////////////////////////////////

/**
 * Frenet Base, defined w.r.t. to a cubic spline
 */
struct FrenetBase2D {
  // Data
  IdType id = InvalidId;  //!< Spline id
  RealType arc_length;    //!< Arc length from start of reference line [m]
  RealType orientation;   //!< Signed orientation of spline at frenet base [rad]
  RealType curvature;     //!< Signed curvature of spline at frenet base [1/m]
  RealType curvature_change_rate;  //!< Signed curvature change rate [1/m^2]
  cubic_spline::SegmentInfo<IdxType, RealType>
      segment_info;  //!< point on segments

  FrenetBase2D()
      : arc_length(-1.0),
        orientation(0.0),
        curvature(0.0),
        curvature_change_rate(0.0),
        segment_info(0, 0.0) {}
  FrenetBase2D(const IdType _id, const RealType _arclength,
               const RealType _orientation, const RealType _curvature,
               const RealType _curvature_change_rate,
               const cubic_spline::SegmentInfo<IdxType, RealType> _segment_info)
      : id(_id),
        arc_length(_arclength),
        orientation(_orientation),
        curvature(_curvature),
        curvature_change_rate(_curvature_change_rate),
        segment_info(_segment_info) {}
  explicit FrenetBase2D(const FrenetBase2D& other)
      : id(other.id),
        arc_length(other.arc_length),
        orientation(other.orientation),
        curvature(other.curvature),
        curvature_change_rate(other.curvature_change_rate),
        segment_info(other.segment_info) {}
};

inline std::ostream& operator<<(std::ostream& os, const FrenetBase2D& fb) {
  os << "Frenet base: " << fb.id << ", l = " << fb.arc_length
     << ", orientation = " << fb.orientation << ", curvature = " << fb.curvature
     << ", CCR = " << fb.curvature_change_rate << " [ " << fb.segment_info.idx
     << ", " << fb.segment_info.relative_arc_length << " ]\n";
  return os;
};

// Forward declaration
class FrenetState2D;
struct FrenetStateVector2D;
struct FrenetStateCovarianceMatrix2D;

/**
 * @brief Frenet Frame, which can be specified for any point on a cubic spline
 *
 */
class FrenetFrame2D {
 public:
  FrenetFrame2D()
      : frenet_base_(),
        origin_(0.0, 0.0),
        tangent_(0.0, 0.0),
        normal_(0.0, 0.0) {
    rotMat_F2C_ = RotationMatrix::Zero();
    rotMat_C2F_ = RotationMatrix::Zero();
  }
  FrenetFrame2D(const FrenetBase2D& other_base, const CartesianPoint2D& origin,
                const CartesianVector2D& tangent,
                const CartesianVector2D& normal)
      : frenet_base_(other_base),
        origin_(origin),
        tangent_(tangent),
        normal_(normal) {
    // Set rotation matrices (only for internal use)
    rotMat_F2C_ << tangent_.x(), normal_.x(), tangent_.y(), normal_.y();
    rotMat_C2F_ << tangent_.x(), tangent_.y(), normal_.x(), normal_.y();
  }

  /**
   * @brief Conversion from Frenet-based point to Cartesian point.
   *        cartesian_vector = RotMat_F2C * relative_vector
   * @param frenet_position: position with respect to the reference line
   * @return CartesianPoint2Ds
   */
  CartesianPoint2D FromFrenetPoint(const FrenetPoint2D& frenet_point) const;

  CartesianVector2D FromFrenetVector(const FrenetVector2D& frenet_vector) const;

  /**
   * @brief Conversion from a Cartesian-based point to the Frenet Frame
   *        r_cart = RotMat_C2F * local_cartesian_point + frenet_base
   * @param cartesian_point
   * @return FrenetPoint2D
   */
  FrenetPoint2D FromCartesianPoint(
      const CartesianPoint2D& cartesian_point) const;

  FrenetVector2D FromCartesianVector(
      const CartesianVector2D& cartesian_vector) const;

  RealType FromCartesianOrientation(const RealType cartesian_orientation) const;

  FrenetStateVector2D FromCartesianStateVector(
      const CartesianStateVector2D& state_vector,
      const bool moving_frenet_frame = false) const;

  FrenetStateCovarianceMatrix2D FromCartesianStateCovarianceMatrix(
      const CartesianStateCovarianceMatrix2D& state_vector_covariance_matrix,
      const bool moving_frenet_frame = false) const;

  FrenetState2D FromCartesianStateTaylorExpansion(
      const CartesianState2D& cartesian_state) const;

  FrenetState2D FromCartesianState(
      const CartesianState2D& cartesian_state,
      const bool moving_frenet_frame = false) const;

  const FrenetBase2D& frenet_base() const { return frenet_base_; }
  const CartesianPoint2D& origin() const { return origin_; }
  const CartesianVector2D& tangent() const { return tangent_; }
  const CartesianVector2D& normal() const { return normal_; }

  const RealType arc_length() const { return frenet_base_.arc_length; }
  const RealType curvature() const { return frenet_base_.curvature; }

  void setId(const IdType id) { frenet_base_.id = id; }

  // Introspection
  friend std::ostream& operator<<(std::ostream& os, const FrenetFrame2D& ff);

 private:
  FrenetBase2D frenet_base_;
  CartesianPoint2D origin_;
  CartesianVector2D tangent_;
  CartesianVector2D normal_;

  using RotationMatrix = Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>;
  RotationMatrix rotMat_F2C_;  //< From Frenet to Cartesian
  RotationMatrix rotMat_C2F_;  //< From Cartesian to Frenet
};
using FrenetFrames2D = std::vector<FrenetFrame2D>;

inline std::ostream& operator<<(std::ostream& os, const FrenetFrame2D& ff) {
  os << ff.frenet_base();
  os << "Origin: " << ff.origin().transpose() << "\n";
  os << "Tangent: " << ff.tangent().transpose() << "\n";
  os << "Normal: " << ff.normal().transpose() << "\n";
  os << "rotMat_C2F:\n" << ff.rotMat_C2F_ << "\n";
  os << "rotMat_F2C:\n" << ff.rotMat_F2C_ << "\n";
  return os;
};

struct FrenetPositionWithFrame {
  FrenetPositionWithFrame(const FrenetFrame2D _frame,
                          const CartesianPoint2D& c_point)
      : frame(_frame) {
    position = frame.FromCartesianPoint(c_point);
  }

  FrenetFrame2D frame;
  FrenetPoint2D position;
};
using FrenetPositionsWithFrames = std::vector<FrenetPositionWithFrame>;

/**
 * @brief Frenet polyline contains (l,d) coordinates, relative to some
 * reference curve specified by its id. (l coordinate is the arc-length to
 * the point's projection along some reference curve, d coordinate is the
 * signed normal distance to the curve at that point).
 *
 */
class FrenetPolyline {
 public:
  /** Data types of each data point (column) in the augmented data matrix */
  enum DataType { kArclength = 0, kDeviation, kSize };

  /** Data matrix of the polyline information. */
  using DataMatrix = Eigen::Matrix<RealType, DataType::kSize, Eigen::Dynamic>;
  using DataRow = Eigen::Matrix<RealType, 1, Eigen::Dynamic>;
  using DataPoint = FrenetPoint2D;

  // // Constr
  FrenetPolyline(const int size = 0) { data_.resize(DataType::kSize, size); }

  /** Number of points in the polyline definition. */
  int size() const { return data_.cols(); };

  /** The point at index idx in the polyline definition. */
  FrenetPoint2D GetPoint(const int idx) const { return data_.col(idx); }
  void SetPoint(const int idx, const DataPoint& pt) { data_.col(idx) = pt; }

  const FrenetPoint2D operator[](const int idx) const { return data_.col(idx); }

  /**
   * @brief This function returns the deviation 'd' coordinate of the (l,d)
   * Frenet coordinates defining the Frenet polyline. Between two nodes linearly
   * interpolated between two closest Frenet polyline points. If the 'l' value
   * is negative or larger than the arc length, then respectively the 'd'
   * coordinate corresponding to the beginning or end of the curve is returned.
   *
   * @param[in] l: l value to query
   * @return 'd' value of the (l,d) Frenet coordinates corresponding to the
   * query 'l' value.
   */
  RealType deviationAt(const RealType query_l) const;

  // Introspection
  friend std::ostream& operator<<(std::ostream& os, const FrenetPolyline& fpl);

 private:
  /**
   * Data Matrix which describes the polyline
   * Rows:    (l,d) coordinates of polyline points, 1,2,3,...,i,...,M (M =
   * number of points) Columns: one sample point m_i = [s_i, d_i]
   */
  DataMatrix data_;
};

inline std::ostream& operator<<(std::ostream& os, const FrenetPolyline& fpl) {
  using namespace std;
  os << "Frenet Polyline: \n";
  os << "IDX\tPT_S\tPT_D\n";
  for (int idx = 0; idx < fpl.data_.cols(); idx++) {
    os << idx << "\t" << std::setprecision(9)
       << fpl.data_(FrenetPolyline::DataType::kArclength, idx) << "\t"
       << fpl.data_(FrenetPolyline::DataType::kDeviation, idx) << std::endl;
  }
  return os;
};

// /////////////////////////////////////////////////////////////////////////////
// Frenet State and covariances
// /////////////////////////////////////////////////////////////////////////////

struct FrenetStateVector2D
    : public Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign> {
  FrenetStateVector2D(void)
      : Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>() {}
  FrenetStateVector2D(const FrenetPoint2D& position,
                      const FrenetVector2D& velocity) {
    (*this) << position.l(), position.d(), velocity.l(), velocity.d();
  }
  FrenetStateVector2D(const RealType pos_l, const RealType pos_d,
                      const RealType vel_l, const RealType vel_d) {
    (*this) << pos_l, pos_d, vel_l, vel_d;
  }

  using Base = Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>;

  // This constructor allows you to construct FrenetStateVector2D from Eigen
  // expressions
  template <typename OtherDerived>
  FrenetStateVector2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 4, 1, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to FrenetStateVector2D
  template <typename OtherDerived>
  FrenetStateVector2D& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Views
  // const Eigen::Ref<const FrenetPoint2D> position() const {
  //   return static_cast<FrenetPoint2D>(this->head<2>());
  // }

  // Copies
  // TODO: make non-mutable view to the block operation
  FrenetPoint2D position() const { return this->head<2>(); }
  FrenetVector2D velocity() const { return this->tail<2>(); }

  // Non-mutable views
  const RealType l() const { return (*this)[0]; }
  const RealType d() const { return (*this)[1]; }
  const RealType vl() const { return (*this)[2]; }
  const RealType vd() const { return (*this)[3]; }

  // Mutable views
  RealType& l() { return (*this)[0]; }
  RealType& d() { return (*this)[1]; }
  RealType& vl() { return (*this)[2]; }
  RealType& vd() { return (*this)[3]; }

  // TODO: make a mutable view on the eigen matrix possible
  // Eigen::Ref<CartesianVector2D> velocity() { return this->segment<2>(2); }
};

struct FrenetCovarianceMatrix2D
    : public Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign> {
  FrenetCovarianceMatrix2D(void)
      : Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>() {}
  FrenetCovarianceMatrix2D(const RealType ll, const RealType dd,
                           const RealType ld = 0.0) {
    (*this) << ll, ld, ld, dd;
  }

  using Base = Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>;

  // This constructor allows you to construct FrenetCovarianceMatrix2D from
  // Eigen expressions
  template <typename OtherDerived>
  FrenetCovarianceMatrix2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 2, 2, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to
  // FrenetCovarianceMatrix2D
  template <typename OtherDerived>
  FrenetCovarianceMatrix2D& operator=(
      const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Non-mutable views
  const RealType ll() const { return (*this)(0, 0); }
  const RealType dd() const { return (*this)(1, 1); }
  const RealType ld() const { return (*this)(0, 1); }
};

struct FrenetStateCovarianceMatrix2D
    : public Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign> {
  FrenetStateCovarianceMatrix2D(void)
      : Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign>() {}
  FrenetStateCovarianceMatrix2D(const FrenetCovarianceMatrix2D& pos,
                                const FrenetCovarianceMatrix2D& vel,
                                const FrenetCovarianceMatrix2D& pos_vel =
                                    FrenetCovarianceMatrix2D::Zero()) {
    this->block<2, 2>(0, 0) = pos;
    this->block<2, 2>(2, 2) = vel;
    this->block<2, 2>(0, 2) = pos_vel;
    this->block<2, 2>(2, 0) = pos_vel.transpose();
  }
  FrenetStateCovarianceMatrix2D(const RealType ll, const RealType dd,
                                const RealType vlvl, const RealType vdvd,
                                const RealType ld = 0.0,
                                const RealType vlvd = 0.0,
                                const RealType lvl = 0.0,
                                const RealType lvd = 0.0,
                                const RealType dvd = 0.0)
      : FrenetStateCovarianceMatrix2D(
            FrenetCovarianceMatrix2D(ll, dd, ld),
            FrenetCovarianceMatrix2D(vlvl, vdvd, vlvd),
            FrenetCovarianceMatrix2D(lvl, dvd, lvd)) {}

  using Base = Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign>;

  // This constructor allows you to construct
  // FrenetStateCovarianceMatrix2D from Eigen expressions
  template <typename OtherDerived>
  FrenetStateCovarianceMatrix2D(const Eigen::MatrixBase<OtherDerived>& other)
      : Eigen::Matrix<RealType, 4, 4, Eigen::DontAlign>(other) {}

  // This method allows you to assign Eigen expressions to
  // FrenetStateCovarianceMatrix2D
  template <typename OtherDerived>
  FrenetStateCovarianceMatrix2D& operator=(
      const Eigen::MatrixBase<OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  // Copies
  FrenetCovarianceMatrix2D position() const { return this->block<2, 2>(0, 0); }
  FrenetCovarianceMatrix2D velocity() const { return this->block<2, 2>(2, 2); }
  FrenetCovarianceMatrix2D pos_vel() const { return this->block<2, 2>(0, 2); }

  // Non-mutable views
  const RealType ll() const { return (*this)(0, 0); }
  const RealType dd() const { return (*this)(1, 1); }
  const RealType ld() const { return (*this)(0, 1); }

  const RealType vlvl() const { return (*this)(2, 2); }
  const RealType vdvd() const { return (*this)(3, 3); }
  const RealType vlvd() const { return (*this)(2, 3); }
};
// /////////////////////////////////////////////////////////////////////////////
// Frenet state (mean and covariance matrix)
// /////////////////////////////////////////////////////////////////////////////

class FrenetState2D {
 public:
  FrenetState2D(void) : mean_(), cov_mat_() {}
  FrenetState2D(
      const FrenetPoint2D& position, const FrenetVector2D& velocity,
      const FrenetCovarianceMatrix2D& cm_position = FrenetCovarianceMatrix2D(),
      const FrenetCovarianceMatrix2D& cm_velocity = FrenetCovarianceMatrix2D(),
      const FrenetCovarianceMatrix2D& cm_pos_vel =
          FrenetCovarianceMatrix2D::Zero())
      : mean_(position, velocity),
        cov_mat_(cm_position, cm_velocity, cm_pos_vel),
        polar_velocity_state_() {}
  FrenetState2D(const FrenetStateVector2D& mean,
                const FrenetStateCovarianceMatrix2D& cov_mat)
      : mean_(mean), cov_mat_(cov_mat), polar_velocity_state_() {}

  // Simple getter
  RealType l() const { return mean_.l(); }
  RealType d() const { return mean_.d(); }
  RealType vl() const { return mean_.vl(); }
  RealType vd() const { return mean_.vd(); }

  FrenetPoint2D position() const { return mean_.position(); }
  FrenetPoint2D velocity() const { return mean_.velocity(); }
  UncertainValue abs_velocity();
  UncertainValue orientation();

  const FrenetStateVector2D& mean() const { return mean_; }
  const FrenetStateCovarianceMatrix2D& covarianceMatrix() const {
    return cov_mat_;
  }

  // Non-const reference
  FrenetStateVector2D& mean() { return mean_; }
  FrenetStateCovarianceMatrix2D& covarianceMatrix() { return cov_mat_; }

  // Introspection
  friend std::ostream& operator<<(std::ostream& os, const FrenetState2D& state);

 private:
  FrenetStateVector2D mean_;
  FrenetStateCovarianceMatrix2D cov_mat_;

  // optional polar interpretation of the velocity vector. Will only be
  // constructed if needed and then cached.
  const PolarStatePtr getPolarVelocityStatePtr();
  PolarStatePtr polar_velocity_state_;
};

// Introspection
inline std::ostream& operator<<(std::ostream& os, const FrenetState2D& state) {
  using namespace std;
  os << "Frenet State Vector: ";
  os << state.mean_.transpose() << "\n";
  os << "Frenet State CovMat:\n";
  os << state.cov_mat_ << "\n";
  if (state.polar_velocity_state_ != nullptr) {
    os << state.polar_velocity_state_;
  }
  return os;
};

}  // namespace corridor