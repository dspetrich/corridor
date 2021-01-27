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
features.corridor_length = 20
features.l = 0  # <-- X
features.sigma_l = 0.5  # <-- variable
features.obj_length_ratio = 0.2  # <-- variable

# Details
l_o = features.corridor_length * features.obj_length_ratio
x_min = -2*l_o
l_c = features.corridor_length
x_max = l_c - x_min

# problem parameters
n_lengths = 1000
n_sigma = 1000

l_range = np.linspace(x_min, x_max, n_lengths,)
sigma_range = np.linspace(0.01, 10, n_sigma,)

xx = np.zeros((n_lengths, n_sigma), dtype='d')
yy = np.zeros((n_lengths, n_sigma), dtype='d')
zz = np.zeros((n_lengths, n_sigma), dtype='d')

# Details
# populate x,y,z arrays
for i in range(n_lengths):
    for j in range(n_sigma):
        xx[i, j] = l_range[i]
        yy[i, j] = sigma_range[j]
        features.l = l_range[i]
        features.sigma_l = sigma_range[j]
        zz[i, j] = corridor.LongitudinalAssignmentConfidence(features)

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
