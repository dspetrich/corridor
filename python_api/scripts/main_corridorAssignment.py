import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats as stats
import math
import corridor


def plot_assignment_function(w_c, w_o, x_max):
    l_x = ()
    l_y = ()
    if w_o == 0:
        l1_x = (-x_max, -0.5*w_c)
        l1_y = (0.0, 0.0)
        l2_x = np.linspace(-0.5*w_c, 0.5*w_c, 2)
        l2_y = 0*l2_x + 1.0
        l3_x = (0.5*w_c, x_max)
        l3_y = (0.0, 0.0)
        l_x = np.concatenate([l1_x, l2_x, l3_x])
        l_y = np.concatenate([l1_y, l2_y, l3_y])
    else:
        m = 1.0/w_o
        b = 0.5*(1+w_c/w_o)
        x_max = (w_c + w_o)
        x_min = -x_max
        if w_c <= w_o:
            l1_x = np.linspace(-0.5*(w_c+w_o), 0, 2)
            l1_y = m*l1_x + b
            l2_x = np.linspace(0.0, 0.5*(w_c+w_o), 2)
            l2_y = -m*l2_x + b
            l_x = np.concatenate([l1_x, l2_x])
            l_y = np.concatenate([l1_y, l2_y])

        else:
            l1_x = np.linspace(-0.5*(w_c+w_o), -0.5*(w_c-w_o), 2)
            l1_y = m*l1_x + b
            l2_x = np.linspace(-0.5*(w_c-w_o), 0.5*(w_c-w_o), 2)
            l2_y = 0*l2_x + 1.0
            l3_x = np.linspace(0.5*(w_c-w_o), 0.5*(w_c+w_o), 2)
            l3_y = -m*l3_x + b
            l_x = np.concatenate([l1_x, l2_x, l3_x])
            l_y = np.concatenate([l1_y, l2_y, l3_y])
    return l_x, l_y


# Figure
# fig = plt.figure()
# gs = gridspec.GridSpec(2, 2)
# ax1 = fig.add_subplot(gs[0, 0])
# ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
# ax3 = fig.add_subplot(gs[1, 0])
# ax4 = fig.add_subplot(gs[1, 1], sharey=ax3)


# 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# Define plain features of the object
features = corridor.CorridorAssignmentFeature()
features.corridor_width = 4
features.d = 0  # <-- X
features.sigma_d = 0.5  # <-- variable
features.obj_width_ratio = 0.5  # <-- variable

w_c = features.corridor_width
x_max = 2*w_c

# problem parameters
nd = 1000
ns = 1000

d = np.linspace(-x_max, x_max, nd,)
sigma_d = np.linspace(0.01, 10, ns,)

xx = np.zeros((nd, ns), dtype='d')
yy = np.zeros((nd, ns), dtype='d')
zz = np.zeros((nd, ns), dtype='d')

# Details
# populate x,y,z arrays
for i in range(nd):
    for j in range(ns):
        xx[i, j] = d[i]
        yy[i, j] = sigma_d[j]
        features.d = d[i]
        features.sigma_d = sigma_d[j]
        zz[i, j] = corridor.LateralConfidence(features)

surf = ax.plot_surface(xx, yy, zz, cmap=cm.jet,
                       linewidth=0, antialiased=False)

# c = ax.pcolormesh(xx, yy, zz, cmap='RdBu')
plt.show()
