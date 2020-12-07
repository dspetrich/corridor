import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats as stats
import math
import corridor


def create_corridor_boundary_pts(w_c, x_max):
    c_x = ()
    c_y = ()
    c1_x = (-x_max, -0.5*w_c)
    c1_y = (0.0, 0.0)
    c2_x = (-0.5*w_c, 0.5*w_c)
    c2_y = (1.0, 1.0)
    c3_x = (0.5*w_c, x_max)
    c3_y = (0.0, 0.0)
    c_x = np.concatenate([c1_x, c2_x, c3_x])
    c_y = np.concatenate([c1_y, c2_y, c3_y])
    return c_x, c_y


def create_inverted_corridor_boundary_pts(w_c, x_max):
    c_x = ()
    c_y = ()
    c1_x = (-x_max, -0.5*w_c)
    c1_y = (2.0, 2.0)
    c2_x = (-0.5*w_c, 0.5*w_c)
    c2_y = (0.0, 0.0)
    c3_x = (0.5*w_c, x_max)
    c3_y = (2.0, 2.0)
    c_x = np.concatenate([c1_x, c2_x, c3_x])
    c_y = np.concatenate([c1_y, c2_y, c3_y])
    return c_x, c_y


def create_assignment_function(w_c, w_o, x_max):
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
fig = plt.figure()
gs = gridspec.GridSpec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1], sharey=ax3)


# ##############################################################################
# Define plain features of the object
# ##############################################################################
features = corridor.CorridorAssignmentFeature()
features.corridor_width = 4
features.d = 0  # <-- X
features.sigma_d = 0.5  # <-- variable
features.obj_width_ratio = 0.5

# Details
w_c = features.corridor_width
x_max = 2*w_c
x_range = np.linspace(-x_max, x_max, 500)

# Draw the corridor boundaries in each supplot
cx, cy = create_inverted_corridor_boundary_pts(w_c, x_max)
ax1.plot(cx, cy, color='gray', linestyle='--', alpha=0.1)
ax1.fill_between(cx, cy, color='gray', alpha=0.1)
ax2.plot(cx, cy, color='gray', linestyle='--', alpha=0.1)
ax2.fill_between(cx, cy, color='gray', alpha=0.1)
ax3.plot(cx, cy, color='gray', linestyle='--', alpha=0.1)
ax3.fill_between(cx, cy, color='gray', alpha=0.1)
ax4.plot(cx, cy, color='gray', linestyle='--', alpha=0.1)
ax4.fill_between(cx, cy, color='gray', alpha=0.1)

# Assignment function
lx, ly = create_assignment_function(w_c, w_c*features.obj_width_ratio, x_max)
ax1.plot(lx, ly, color='black')
ax1.fill_between(lx, ly, color='black', alpha=0.1)


# ##############################################################################
# Fixed w_c, variable sigma
# ##############################################################################
sigma = np.linspace(0.01, 2, 10)

for s in sigma:
    features.sigma_d = s
    # sample
    mu = 0.0
    sigma = features.sigma_d
    x = np.linspace(mu - 6*sigma, mu + 6*sigma, 100)
    ax1.plot(x, stats.norm.pdf(x, mu, sigma))

    lat_conf = []
    for x in x_range:
        features.d = x
        # for s in sigma:
        lat_conf.append(corridor.LateralAssignmentConfidence(features))
    ax2.plot(x_range, lat_conf)

# ##############################################################################
# Fixed sigma, variable w_c
# ##############################################################################

w_o_ratio_1 = np.linspace(0.01, 5, 20)
# w_o_ratio_2 = np.linspace(1, 1, 10)
w_o_ratio = np.concatenate([w_o_ratio_1])
features.sigma_d = 0.5  # <-- variable

mu = 0.0
sigma = features.sigma_d
x = np.linspace(mu - 6*sigma, mu + 6*sigma, 100)
ax3.plot(x, stats.norm.pdf(x, mu, sigma), color='black')
ax3.fill_between(x, stats.norm.pdf(x, mu, sigma), color='black', alpha=0.1)

for w in w_o_ratio:
    features.obj_width_ratio = w
    lat_conf = []
    for x in x_range:
        features.d = x
        # for s in sigma:
        lat_conf.append(corridor.LateralAssignmentConfidence(features))
    ax4.plot(x_range, lat_conf)
    lx, ly = create_assignment_function(w_c, w_c*w, x_max)
    ax3.plot(lx, ly)

# ##############################################################################
# Plot setup
# ##############################################################################
ax1.set_ylim([0, 1.1])
ax1.set_xlim([-x_max, x_max])
ax2.set_xlim([-x_max, x_max])
ax3.set_ylim([0, 1.1])
ax3.set_xlim([-x_max, x_max])
ax4.set_xlim([-x_max, x_max])

plt.setp(ax2.get_yticklabels(), visible=True)
plt.setp(ax4.get_yticklabels(), visible=True)

plt.setp([ax1], title='Fixed object width, variable standard deviation')
plt.setp([ax3], title='Fixed standard deviation, variable object width')
plt.setp([ax2], title='Assignment likelihood function')
plt.setp([ax4], title='Assignment likelihood function')
fig.suptitle('Lateral Assignment', size=20)
gs.tight_layout(fig, rect=[0, 0, 1, 0.97])

plt.show()
