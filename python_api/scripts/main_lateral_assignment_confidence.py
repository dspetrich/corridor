import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
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

# 3D plot
fig = plt.figure(figsize=plt.figaspect(0.5))
# set up the axes for the first plot
ax_sigma_d = fig.add_subplot(1, 2, 1, projection='3d')
ax_width = fig.add_subplot(1, 2, 2, projection='3d')


# Define plain features of the object
features = corridor.CorridorAssignmentFeature()
features.corridor_width = 4
features.d = 0  # <-- X
features.sigma_d = 0.5  # <-- variable
features.obj_width_ratio = 0.5  # <-- variable

w_c = features.corridor_width
x_max = 2.5*w_c

# problem parameters
n_d = 1000
n_sd = 1000
n_w = 1000

d = np.linspace(-x_max, x_max, n_d,)
sigma_d = np.linspace(0.01, 10, n_sd,)
w = np.linspace(0.0, 3, n_w,)


# variable sigma d
xx_sd = np.zeros((n_d, n_sd), dtype='d')
yy_sd = np.zeros((n_d, n_sd), dtype='d')
zz_sd = np.zeros((n_d, n_sd), dtype='d')

# variable relative object width
xx_w = np.zeros((n_d, n_w), dtype='d')
yy_w = np.zeros((n_d, n_w), dtype='d')
zz_w = np.zeros((n_d, n_w), dtype='d')

# Details
# populate x,y,z arrays
for i in range(n_d):
    # variable sigma_d
    features.obj_width_ratio = 0.5
    for j in range(n_sd):
        xx_sd[i, j] = d[i]
        yy_sd[i, j] = sigma_d[j]
        features.d = d[i]
        features.sigma_d = sigma_d[j]
        zz_sd[i, j] = corridor.LateralAssignmentConfidence(features)
    # variable width
    features.sigma_d = 0.5
    for j in range(n_w):
        xx_w[i, j] = d[i]
        yy_w[i, j] = w[j]
        features.d = d[i]
        features.obj_width_ratio = w[j]
        zz_w[i, j] = corridor.LateralAssignmentConfidence(features)

surf_sd = ax_sigma_d.plot_surface(xx_sd, yy_sd, zz_sd, cmap=cm.coolwarm,
                                  linewidth=0, antialiased=False)
surf_w = ax_width.plot_surface(xx_w, yy_w, zz_w, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

ax_sigma_d.set_xlabel('lateral position $d_{\eta}$ [m]')
ax_sigma_d.set_ylabel('standard deviation $\sigma_{d}$ [m]')
ax_sigma_d.set_zlabel('assignment confidence')

ax_width.set_xlabel('lateral position $d_{\eta}$ [m]')
ax_width.set_ylabel('object width ratio $\hat{W}_{obj}/W_{corr}$')
ax_width.set_zlabel('assignment confidence')

# cset = ax.contour(xx_sd, yy_sd, zz_sd, zdir='x', cmap=cm.coolwarm)
# cset = ax.contour(xx_sd, yy_sd, zz_sd, zdir='y', cmap=cm.coolwarm)

# c = ax.pcolormesh(xx_sd, yy_sd, zz_sd, cmap='RdBu')
# plt.savefig(
#     '/home/dsp/Pictures/Matplotlib_PGFs/CorridorAssignment.pdf', bbox_inches='tight')
plt.savefig(
    '/home/dsp/Pictures/Matplotlib_PGFs/LateralAssignmentConfidence.pdf')
plt.show()
