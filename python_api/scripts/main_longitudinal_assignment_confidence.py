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
ax_sigma_l = fig.add_subplot(1, 2, 1, projection='3d')
ax_length = fig.add_subplot(1, 2, 2, projection='3d')


# Define plain features of the object
features = corridor.CorridorAssignmentFeature()
features.corridor_length = 10
features.l = 0  # <-- X
features.sigma_l = 0.5  # <-- variable
features.obj_length_ratio = 0.2  # <-- variable

# Details
l_obj = features.corridor_length * features.obj_length_ratio
x_min = -2*l_obj
l_c = features.corridor_length
x_max = l_c - x_min

# problem parameters
n_l = 1000
n_sd = 1000
n_w = 1000

l = np.linspace(x_min, x_max, n_l,)
sigma_l = np.linspace(0.01, 10, n_sd,)
w = np.linspace(0.0, 2, n_w,)

# variable sigma d
xx_sd = np.zeros((n_l, n_sd), dtype='d')
yy_sd = np.zeros((n_l, n_sd), dtype='d')
zz_sd = np.zeros((n_l, n_sd), dtype='d')

# variable relative object width
xx_w = np.zeros((n_l, n_w), dtype='d')
yy_w = np.zeros((n_l, n_w), dtype='d')
zz_w = np.zeros((n_l, n_w), dtype='d')

# Details
# populate x,y,z arrays
for i in range(n_l):
    # variable sigma_l
    features.obj_length_ratio = 0.2
    for j in range(n_sd):
        xx_sd[i, j] = l[i]
        yy_sd[i, j] = sigma_l[j]
        features.l = l[i]
        features.sigma_l = sigma_l[j]
        zz_sd[i, j] = corridor.LongitudinalAssignmentConfidence(features)
    # variable width
    features.sigma_l = 0.5
    for j in range(n_w):
        xx_w[i, j] = l[i]
        yy_w[i, j] = w[j]
        features.d = l[i]
        features.obj_length_ratio = w[j]
        zz_w[i, j] = corridor.LongitudinalAssignmentConfidence(features)

surf_sd = ax_sigma_l.plot_surface(xx_sd, yy_sd, zz_sd, cmap=cm.coolwarm,
                                  linewidth=0, antialiased=False)
surf_w = ax_length.plot_surface(xx_w, yy_w, zz_w, cmap=cm.coolwarm,
                                linewidth=0, antialiased=False)

ax_sigma_l.set_xlabel('longitudinal position $l_{r}$ [m]')
ax_sigma_l.set_ylabel('standard deviation $\sigma_{l}$ [m]')
ax_sigma_l.set_zlabel('assignment confidence')

ax_length.set_xlabel('longitudinal position $l_{r}$ [m]')
ax_length.set_ylabel('object length ratio $\hat{l}_{obj}/l_{corr}$')
ax_length.set_ylabel('object length ratio $\hat{l}_{obj}/l_{corr}$')
ax_length.set_zlabel('assignment confidence')

# cset = ax.contour(xx_sd, yy_sd, zz_sd, zdir='x', cmap=cm.coolwarm)
# cset = ax.contour(xx_sd, yy_sd, zz_sd, zdir='y', cmap=cm.coolwarm)

# c = ax.pcolormesh(xx_sd, yy_sd, zz_sd, cmap='RdBu')
# plt.savefig(
#     '/home/dsp/Pictures/Matplotlib_PGFs/CorridorAssignment.pdf', bbox_inches='tight')
plt.savefig(
    '/tmp/LongitudinalAssignmentConfidence.pdf')
plt.show()
