import random
import statistics
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import corridor


def sample_polar_data(mean, covMat, n_samples=5000):
    r, phi = np.random.multivariate_normal(mean, covMat, n_samples).T
    r = abs(r)
    np.unwrap(phi)
    return np.stack((r, phi), axis=0)


def confidence_ellipse(ax, mean, cov_mat, n_std=3, facecolor='none', **kwargs):
    pearson = cov_mat[0, 1]/np.sqrt(cov_mat[0, 0] * cov_mat[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov_mat[0, 0]) * n_std
    mean_x = mean[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov_mat[1, 1]) * n_std
    mean_y = mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


# Standard deviation for error ellipses
n_std = 2

# Figure
fig, (ax_init_polar, ax_cart, ax_polar) = plt.subplots(1, 3)

# Initial velocity value and orientation
radius = 10
heading = 1 * math.pi / 4.0

var_radius = 0.5
var_heading = 1 * math.pi / 32.0
cov_radHead = 0

# Initial 2d p vector with normal distribution
initial_polar_mean = [radius, heading]
initial_polar_covMat = [
    [var_radius, cov_radHead], [cov_radHead, var_heading]]

# Convert to cartesian and calculate new mean and cov
initial_polar_data = sample_polar_data(
    initial_polar_mean, initial_polar_covMat, n_samples=5000)

# mean and covMat from data
polar_mean_0 = np.mean(initial_polar_data, axis=1)
polar_cov_0 = np.cov(initial_polar_data)

# Initial subplot
ax_init_polar.plot(
    initial_polar_data[0, :], initial_polar_data[1, :], 'g.', alpha=0.3)
ax_init_polar.plot(polar_mean_0[0], polar_mean_0[1], 'ro')
confidence_ellipse(ax_init_polar, polar_mean_0, polar_cov_0,
                   n_std, facecolor='none', zorder=10, edgecolor='r')

# Monte Carlo Methode for polar to cartesian transformation
x_list = []
y_list = []
for i in range(np.size(initial_polar_data, 1)):
    result = corridor.polar_to_cartesian_2d(
        initial_polar_data[0, i], initial_polar_data[1, i])
    x_list.append(result[0])
    y_list.append(result[1])

# mean and covMat from data
cart_data = np.stack((x_list, y_list), axis=0)
cart_mean = np.mean(cart_data, axis=1)
cart_cov = np.cov(cart_data)

ax_cart.plot(x_list, y_list, 'g.', alpha=0.3)
ax_cart.plot(cart_mean[0], cart_mean[1], 'ro')
confidence_ellipse(ax_cart, cart_mean, cart_cov,
                   n_std, facecolor='none', zorder=10, edgecolor='r')

# Monte Carlo Methode for cartesian to polar transformation
r_list = []
phi_list = []
for i in range(np.size(cart_data, 1)):
    result = corridor.cartesian_to_polar_2d(
        cart_data[0, i], cart_data[1, i])
    r_list.append(result[0])
    phi_list.append(result[1])

# mean and covMat from data
polar_data = np.stack((r_list, phi_list), axis=0)
polar_mean = np.mean(polar_data, axis=1)
polar_cov = np.cov(polar_data)

ax_polar.plot(r_list, phi_list, 'g.', alpha=0.3)
ax_polar.plot(polar_mean[0], polar_mean[1], 'ro')
confidence_ellipse(ax_polar, polar_mean, polar_cov,
                   n_std, facecolor='none', zorder=10, edgecolor='r')

# Unscented tranformation from cartesian to polar
cart_state = corridor.FlatCartesianStateAndCovMat2D()
cart_state.x = cart_mean[0]
cart_state.y = cart_mean[1]
cart_state.var_x = cart_cov[0, 0]
cart_state.var_y = cart_cov[1, 1]
cart_state.cov_xy = cart_cov[1, 0]

polar_state = corridor.ut_cartesian_to_polar_2d(cart_state)

# Create mean and cov mat
ut_polar_mean = np.array([polar_state.r, polar_state.phi])
ut_polar_cov = np.array([[polar_state.var_r, polar_state.cov_rphi], [
    polar_state.cov_rphi, polar_state.var_phi]])

ax_polar.plot(ut_polar_mean[0], ut_polar_mean[1], 'bo')
confidence_ellipse(ax_polar, ut_polar_mean, ut_polar_cov,
                   n_std, facecolor='none', zorder=10, edgecolor='b', linewidth=2)

ax_cart.axis('equal')

plt.show()
