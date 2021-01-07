# Check UT transformation for polar coordinate tranformation against MC

import random
import statistics
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import corridor

matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': '10',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.autolayout': True,
    # 'figure.figsize': [7, 4],
    'axes.titlesize': 'medium',
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    'legend.fontsize': 'x-small',
    'legend.title_fontsize': 'small',
    # 'axes.labelsize': 'small',
})


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


# Initial 2d cartesian vector with normal distribution
mean = [5, 5]
covMat = [[1, 0.5],
          [0.5, 2]]


# Draw cartesian samples
n_samples = 5000
x, y = np.random.multivariate_normal(mean, covMat, n_samples).T

# Calculate mean and covariance matrix of the sample data
cart_data = np.stack((x, y), axis=0)
cart_mean = np.mean(cart_data, axis=1)
cart_cov = np.cov(cart_data)

# print("Cartesian State")
# print(cart_mean)
# print(cart_cov)

fig, (ax_cart, ax_polar) = plt.subplots(1, 2)

ax_cart.plot(x, y, 'g.', alpha=0.3)
ax_cart.axis('equal')

# Plot mean and covMat
ax_cart.plot(mean[0], mean[1], 'ro')
confidence_ellipse(ax_cart, cart_mean, cart_cov,
                   3, facecolor='none', zorder=10, edgecolor='r')


# Monte Carlo Methode for polar coordinate transformation
r_list = []
phi_list = []
for i in range(n_samples):
    result = corridor.cartesian_to_polar_2d(x[i], y[i])
    r_list.append(result[0])
    phi_list.append(result[1])

ax_polar.plot(r_list, phi_list, 'g.', alpha=0.3)

# Calculate mean and covariance matrix
polar = np.stack((r_list, phi_list), axis=0)
polar_mean = np.mean(polar, axis=1)
polar_cov = np.cov(polar)

# print("Polar State")
# print(polar_mean)
# print(polar_cov)


ax_polar.plot(polar_mean[0], polar_mean[1], 'ro')
confidence_ellipse(ax_polar, polar_mean, polar_cov,
                   3, facecolor='none', zorder=10, edgecolor='r')


# Unscented tranformation
cart_state = corridor.FlatCartesianPositionAndCovMat2D()
cart_state.x = cart_mean[0]
cart_state.y = cart_mean[1]
cart_state.var_x = cart_cov[0, 0]
cart_state.var_y = cart_cov[1, 1]
cart_state.cov_xy = cart_cov[1, 0]

polar_state = corridor.ut_cartesian_to_polar_2d(cart_state)

print("UT Polar State")
print(polar_state.r, polar_state.phi)

# Create mean and cov mat
ut_polar_mean = np.array([polar_state.r, polar_state.phi])
ut_polar_cov = np.array([[polar_state.var_r, polar_state.cov_rphi], [
    polar_state.cov_rphi, polar_state.var_phi]])

# print(ut_polar_mean)
# print(ut_polar_cov)

ax_polar.plot(ut_polar_mean[0], ut_polar_mean[1], 'bo')
confidence_ellipse(ax_polar, ut_polar_mean, ut_polar_cov,
                   3, facecolor='none', zorder=10, edgecolor='b')


# r_mean = sum(r_list) / len(r_list)
# phi_mean = sum(phi_list) / len(phi_list)
# print(r_mean)
# print(phi_mean)
# print(statistics.mean(r_list))
# print(statistics.mean(phi_list))

plt.show()
