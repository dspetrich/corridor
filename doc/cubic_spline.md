# Cubic Spline

The cubic spline implementation can be found [here](../include/corridor/cubic_spline/.).

## Overview
The purpose of the cubic spline generation is to provide a parametric curve description utilized to define a Frenet frame at any arc-length $l$.

At the moment, it is realized as a natural cubic splu

## Spline creation
- from polygon to cubic spline

## Coordinate transformation

- Transformation functions for position and velocity
- Error propagation for the two assumptions of the projection point.
  - Linear projection, when assumed that the projection point is not in motion
  - non-linear projection in case the projection point has the same velocity as the tangential velocity of the source.
  - $E=m c^2 \qquad \text{with} \quad \frac{1}{2}\pi$