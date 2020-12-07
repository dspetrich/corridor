import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats as stats
import math
import corridor


# 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# Define plain features of the object
# problem parameters
n_vel = 1000
n_std = 1000

vel = np.linspace(0, 20, n_vel,)
sigma_vel = np.linspace(1e-6, 10, n_std,)

xx = np.zeros((n_vel, n_std), dtype='d')
yy = np.zeros((n_vel, n_std), dtype='d')
zz = np.zeros((n_vel, n_std), dtype='d')

# Details
# populate x,y,z arrays
for i in range(n_vel):
    for j in range(n_std):
        xx[i, j] = vel[i]
        yy[i, j] = sigma_vel[j]
        zz[i, j] = corridor.MovingConfidence(vel[i], sigma_vel[j], 3.0)

surf = ax.plot_surface(xx, yy, zz, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True,
                       # rstride=10, cstride=10,
                       alpha=0.8)
# cset = ax.contour(xx, yy, zz, zdir='z', offset=1, cmap=cm.coolwarm)
cset = ax.contour(xx, yy, zz, zdir='x', offset=20, cmap=cm.coolwarm)
cset = ax.contour(xx, yy, zz, zdir='y', offset=10, cmap=cm.coolwarm)

ax.set_xlabel('absolute velocity [mps]')
ax.set_ylabel(r'$\sigma_{velocity}$ [mps]')
ax.set_zlabel('Moving Confidence')

plt.show()
