# Corridor

## Overview

Corridor is a C++ library for creating a curve-based description of a traffic lane segment, utilizing the [Frenet–Serret formulas](https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas). Such representations can be applied in the context of automated driving to describe traffic situations within the topology of the road network.
The library is designed as an extension to the [lanelet2](https://github.com/fzi-forschungszentrum-informatik/Lanelet2) map library but could potentially be used with any other map format.
<!-- , as long the required inputs are provided for creating a corridor element. -->

### Features:
- Creates a **corridor from three 2D polygons**, representing the left and right lane boundary, and the reference line in a Cartesian coordinate frame (e.g. UTM or similar).
- Reference line is converted into a **2D cubic spline**, which serves as reference curve of the Frenet frame. Besides being located between the left and right boundary, no further requirements towards the reference polygon need to be meet (no centerline). However, it is recommended that the sampling rate of the polygon is higher in areas with high curvature to ensure good interpolation behavior. All three polygons can have their own sampling rate which can variate, e.g. no equidistant sampling required. 
- Left and right boundary are directly converted into a Frenet polygon representation, consisting of arc length and deviation with respect to the reference line. This allows an efficient query of the corridor width and signed boundary distances at any query length.
- Consecutive corridors can be combined to a continuos **corridor sequence**, which behaves similar to a single corridor.
- **Utility functions** for coordinate transformations under error propagation for linear and non-linear projection. Build-in projections methods for position and velocity information including the covariance matrix and polar coordinate transformation (particular usefully for velocity information)
- Non-linear projection based on a **generic Unscented Transformation** implementation. Details can be found here:
- **Probabilistic assignment functions** for any objects and points to the corridor. This includes relative location and orientation of the object/point to the corridor. Furter information in the documentation:
- **Python** bindings and evaluation scripts for the most important queries.
- Released under the [**BSD 3-Clause license**](LICENSE)

![](lanelet2_core/doc/images/lanelet2_example_image.png)

## Documentation

You can find more documentation in the documentation and in doxygen comments within the code. Here is an overview on the most important topics:

- [Here](doc/cubic_spline.md) you'll find more information about the cubic spline creation, the underlying math and the coordinate transformation with error propagation.

## Installation

Corridor relies mainly on [Catkin](https://catkin-tools.readthedocs.io/en/latest/index.html) for building and is targeted towards Linux.

At least **C++14** is required.

### Dependencies
Besides [Catkin](https://catkin-tools.readthedocs.io/en/latest/index.html), the dependencies are
* `Boost` (from 1.58)
* `eigen3`
* [`mrt_cmake_modules`](https://github.com/KIT-MRT/mrt_cmake_modules), a CMake helper library
* `boost-python, python2 or python3` (for python_api)

For Ubuntu, the steps are the following:
* [Set up ROS](http://wiki.ros.org/ROS/Installation), and install at least `rospack`, `catkin` and `mrt_cmake_modules` (e.g. `ros-melodic-rospack`, `ros-melodic-catkin`, `ros-melodic-mrt-cmake-modules`):
```
sudo apt-get install ros-melodic-rospack ros-melodic-catkin ros-melodic-mrt-cmake-modules
```

* Install the dependencies above:
```bash
sudo apt-get install libboost-dev libeigen3-dev libpugixml-dev libpython-dev libboost-python-dev python-catkin-tools
```

**On 16.04 and below**, `mrt_cmake_modules` is not available in ROS and you have to clone it into your workspace (`git clone https://github.com/KIT-MRT/mrt_cmake_modules.git`).

### Building
As usual with Catkin, after you have sourced the ros installation, you have to create a workspace and clone all required packages there. Then you can build.
```shell
source /opt/ros/$ROS_DISTRO/setup.bash
mkdir catkin_ws && cd catkin_ws && mkdir src
catkin init
catkin config --cmake-args -DCMAKE_BUILD_TYPE=RelWithDebInfo # build in release mode (or whatever you prefer)
cd src
git clone https://THIS_REPO
cd ..
catkin build
```

### Python3

The python bindings are implemented and tested with python 3. In general it should work with pythion 2.7 and above as well, but I didn't test it yet.

## Citation

If you are using Corridor for scientific research or work, it would be nice if the usage of this library is mentioned. A scientific publication of this work and its usage is currently in the making and will be mentioned here shortly.


