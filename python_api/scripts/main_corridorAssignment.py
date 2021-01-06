import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
        zz[i, j] = corridor.LateralAssignmentConfidence(features)

surf = ax.plot_surface(xx, yy, zz, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlabel('Lateral displacement $d_{p,r}$ [$m$]')
ax.set_ylabel('Standard Deviation $\sigma_{d}$ [$m$]')
ax.set_zlabel('Assignment confidence')

cset = ax.contour(xx, yy, zz, zdir='x', cmap=cm.coolwarm)
cset = ax.contour(xx, yy, zz, zdir='y', cmap=cm.coolwarm)

# c = ax.pcolormesh(xx, yy, zz, cmap='RdBu')
# plt.savefig(
#     '/home/dsp/Pictures/Matplotlib_PGFs/CorridorAssignment.pdf', bbox_inches='tight')
plt.savefig(
    '/home/dsp/Pictures/Matplotlib_PGFs/CorridorAssignment.pdf')
plt.show()
