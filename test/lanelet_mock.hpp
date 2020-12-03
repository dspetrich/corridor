#pragma once

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

struct FirstStraightLaneletSection : public LaneletMock {
  void fillCenterline() {
    centerline.emplace_back(0.f, 0.f);
    centerline.emplace_back(5.f, 0.f);
    centerline.emplace_back(10.f, 0.f);
  }
  void fillLeftBoundary() {
    left_boundary.emplace_back(0, 2.f);
    left_boundary.emplace_back(1.f, 2.f);
    left_boundary.emplace_back(2.f, 2.f);
    left_boundary.emplace_back(3.f, 2.f);
    left_boundary.emplace_back(4.f, 2.f);
    left_boundary.emplace_back(5.f, 2.f);
    left_boundary.emplace_back(6.f, 2.f);
    left_boundary.emplace_back(7.f, 2.f);
    left_boundary.emplace_back(8.f, 2.f);
    left_boundary.emplace_back(9.f, 2.f);
    left_boundary.emplace_back(10.f, 2.f);
  }
  void fillRightBoundary() {
    right_boundary.emplace_back(0, -2.f);
    right_boundary.emplace_back(1.f, -2.f);
    right_boundary.emplace_back(2.f, -2.f);
    right_boundary.emplace_back(3.f, -2.f);
    right_boundary.emplace_back(4.f, -2.f);
    right_boundary.emplace_back(5.f, -2.f);
    right_boundary.emplace_back(6.f, -2.f);
    right_boundary.emplace_back(7.f, -2.f);
    right_boundary.emplace_back(8.f, -2.f);
    right_boundary.emplace_back(9.f, -2.f);
    right_boundary.emplace_back(10.f, -2.f);
  }
  FirstStraightLaneletSection() {
    fillCenterline();
    fillLeftBoundary();
    fillRightBoundary();
  }
};

struct SecondStraightLaneletSection : public LaneletMock {
  void fillCenterline() {
    centerline.emplace_back(10.f, 0.f);
    centerline.emplace_back(15.f, 0.f);
    centerline.emplace_back(20.f, 0.f);
  }
  void fillLeftBoundary() {
    left_boundary.emplace_back(10, 2.f);
    left_boundary.emplace_back(11.f, 2.f);
    left_boundary.emplace_back(12.f, 2.f);
    left_boundary.emplace_back(13.f, 2.f);
    left_boundary.emplace_back(14.f, 2.f);
    left_boundary.emplace_back(15.f, 2.f);
    left_boundary.emplace_back(16.f, 2.f);
    left_boundary.emplace_back(17.f, 2.f);
    left_boundary.emplace_back(18.f, 2.f);
    left_boundary.emplace_back(19.f, 2.f);
    left_boundary.emplace_back(20.f, 2.f);
  }
  void fillRightBoundary() {
    right_boundary.emplace_back(10, -2.f);
    right_boundary.emplace_back(11.f, -2.f);
    right_boundary.emplace_back(12.f, -2.f);
    right_boundary.emplace_back(13.f, -2.f);
    right_boundary.emplace_back(14.f, -2.f);
    right_boundary.emplace_back(15.f, -2.f);
    right_boundary.emplace_back(16.f, -2.f);
    right_boundary.emplace_back(17.f, -2.f);
    right_boundary.emplace_back(18.f, -2.f);
    right_boundary.emplace_back(19.f, -2.f);
    right_boundary.emplace_back(20.f, -2.f);
  }
  SecondStraightLaneletSection() {
    fillCenterline();
    fillLeftBoundary();
    fillRightBoundary();
  }
};