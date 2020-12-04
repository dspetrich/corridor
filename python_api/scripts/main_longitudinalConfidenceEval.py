import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats as stats
import math
import corridor


def create_assignment_function(l_c, l_o, x_min, x_max):
    l_x = ()
    l_y = ()
    if l_o == 0:
        l1_x = (x_min, 0.0)
        l1_y = (0.0, 0.0)
        l2_x = np.linspace(0.0, l_c, 2)
        l2_y = 0*l2_x + 1.0
        l3_x = (l_c, x_max)
        l3_y = (0.0, 0.0)
        l_x = np.concatenate([l1_x, l2_x, l3_x])
        l_y = np.concatenate([l1_y, l2_y, l3_y])
    else:
        m1 = 1.0/l_o
        b1 = 0.5
        m2 = -m1
        b2 = 0.5 + l_c/l_o

        if l_c <= l_o:
            l1_x = np.linspace(-0.5*l_o, 0.5*l_c, 2)
            l1_y = m1*l1_x + b1
            l2_x = np.linspace(0.5*l_c, l_c+0.5*l_o, 2)
            l2_y = m2*l2_x + b2
            l_x = np.concatenate([l1_x, l2_x])
            l_y = np.concatenate([l1_y, l2_y])

        else:
            l1_x = np.linspace(-0.5*l_o, 0.5*l_o, 2)
            l1_y = m1*l1_x + b1
            l2_x = np.linspace(0.5*l_o, l_c-0.5*l_o, 2)
            l2_y = 0*l2_x + 1.0
            l3_x = np.linspace(l_c-0.5*l_o, l_c+0.5*l_o, 2)
            l3_y = m2*l3_x + b2
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

x_range = np.linspace(x_min, x_max, 1000)

# Assignment function
lx, ly = create_assignment_function(l_c, l_o, x_min, x_max)
ax1.plot(lx, ly, color='black')
ax1.fill_between(lx, ly, color='black', alpha=0.1)

# Fixed l_c, variable sigma
sigma = np.linspace(0.01, 5, 10)

for l in sigma:
    features.sigma_l = l
    # sample
    mu = 0.5*l_c
    sigma = features.sigma_l
    x = np.linspace(mu - 6*sigma, mu + 6*sigma, 100)
    ax1.plot(x, stats.norm.pdf(x, mu, sigma))

    long_conf = []
    for x in x_range:
        features.l = x
        # for s in sigma:
        long_conf.append(corridor.LongitudinalConfidence(features))
    ax2.plot(x_range, long_conf)


# Fixed sigma, variable l_c
l_c_ratio = np.linspace(0.001, 1, 20)
features.sigma_l = 0.5  # <-- variable

mu = 0.5*l_c
sigma = features.sigma_l
x = np.linspace(mu - 6*sigma, mu + 6*sigma, 100)
ax3.plot(x, stats.norm.pdf(x, mu, sigma), color='black')
ax3.fill_between(x, stats.norm.pdf(x, mu, sigma), color='black', alpha=0.1)

for l in l_c_ratio:
    features.obj_length_ratio = l
    long_conf = []
    for x in x_range:
        features.l = x
        # for s in sigma:
        long_conf.append(corridor.LongitudinalConfidence(features))
    ax4.plot(x_range, long_conf)
    lx, ly = create_assignment_function(l_c, l_c*l, x_min, x_max)
    ax3.plot(lx, ly)

# Plot setup
ax1.set_ylim([0, 1.1])
ax3.set_ylim([0, 1.1])

ax1.set_xlim([x_min, x_max])
ax2.set_xlim([x_min, x_max])
ax3.set_xlim([x_min, x_max])
ax4.set_xlim([x_min, x_max])

plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)

plt.setp([ax1], title='Fixed object width, variable standard deviation')
plt.setp([ax3], title='Fixed standard deviation, variable object width')
plt.setp([ax2], title='Assignment likelihood function')
plt.setp([ax4], title='Assignment likelihood function')
fig.suptitle('Longitudinal Assignment', size=20)
gs.tight_layout(fig, rect=[0, 0, 1, 0.97])

plt.show()
