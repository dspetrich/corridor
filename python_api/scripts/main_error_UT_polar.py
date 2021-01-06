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


def approximation_error(radius, var_radius, heading, var_heading, cov_radHead):
    # Initial 2d p vector with normal distribution
    initial_polar_mean = [radius, heading]
    initial_polar_covMat = [
        [var_radius, cov_radHead], [cov_radHead, var_heading]]

    # Extract samples from distribution
    initial_polar_data = sample_polar_data(
        initial_polar_mean, initial_polar_covMat, n_samples=5000)

    # mean and covMat from data
    polar_mean_0 = np.mean(initial_polar_data, axis=1)
    polar_cov_0 = np.cov(initial_polar_data)

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

    # Unscented tranformation from cartesian to polar
    cart_state = corridor.FlatCartesianPositionAndCovMat2D()
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

    delta = abs(polar_state.r - polar_mean[0])
    return delta


n_std = 2

fig, ax = plt.subplots()

r_range = np.linspace(0, 20, 40)
std_r_range = np.linspace(0.5, 10, 10)

# Heading angle doesn't play a big role in the error.
heading_range = np.linspace(0.0, 2*math.pi, 4)

# But of cause the standard deviation of the heading
std_h_range = np.linspace(1e-5, math.pi/4.0, 10)

lines = []
for std_h in std_h_range:
    for std_r in std_r_range:
        data = []
        for r in r_range:
            # Initial velocity value and orientation
            radius = r
            heading = math.pi/4.0

            var_radius = std_r * std_r
            var_heading = std_h * std_h
            cov_radHead = 0

            data.append(approximation_error(
                r, var_radius, heading, var_heading, cov_radHead))

        lines.append(data)
        ax.plot(r_range, data, '-', linewidth=1)

plt.show()
