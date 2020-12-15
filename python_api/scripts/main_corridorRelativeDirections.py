import matplotlib
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from matplotlib.ticker import FuncFormatter, MultipleLocator, FormatStrFormatter
import numpy as np

import corridor


# matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': '10',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.autolayout': True,
    'figure.figsize': [5, 3],
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    'legend.fontsize': 'x-small',
    'legend.title_fontsize': 'small',
    'lines.linewidth': '0.5'
    # 'axes.labelsize': 'small',
})

# ##############################################################################
# Axis formatting with fractions of pi (or other constants)
# ##############################################################################
# source: https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlibs


def multiple_formatter(denominator=4, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num, den)
        (num, den) = (int(num/com), int(den/com))
        if den == 1:
            if num == 0:
                return r'$0$'
            if num == 1:
                return r'$%s$' % latex
            elif num == -1:
                return r'$-%s$' % latex
            else:
                return r'$%s%s$' % (num, latex)
        else:
            if num == 1:
                return r'$\frac{%s}{%s}$' % (latex, den)
            elif num == -1:
                return r'$\frac{-%s}{%s}$' % (latex, den)
            else:
                return r'$\frac{%s%s}{%s}$' % (num, latex, den)
    return _multiple_formatter


def num_pi_formatter(x):
    num = x/np.pi
    if num < 1e-3:
        return '{:0.1e}'.format(num) + r' $\pi$'
        return str(num) + r' $\pi$'
    else:
        return str(round(num, 2)) + r' $\pi$'


# ##############################################################################
# Data generation
# ##############################################################################
# Set samples for heading and its standard deviation
n_heading_angles = 1000
n_heading_std = 6

x_max = 1.*np.pi
x_min = -x_max

heading_angle = np.linspace(x_min, x_max, n_heading_angles,)
sigma_heading_angle = np.linspace(1e-12, 0.25, n_heading_std,)*np.pi

downstream_lines = []
upstream_lines = []
towards_left_lines = []
towards_right_lines = []
for i in range(n_heading_std):
    downstream = []
    upstream = []
    towards_left = []
    towards_right = []
    for j in range(n_heading_angles):
        result_dict = corridor.RelativeDirectionConfidence(
            heading_angle[j], sigma_heading_angle[i], 3)
        downstream.append(result_dict['downstream'])
        upstream.append(result_dict['upstream'])
        towards_left.append(result_dict['towards_left'])
        towards_right.append(result_dict['towards_right'])

    downstream_lines.append(downstream)
    upstream_lines.append(upstream)
    towards_left_lines.append(towards_left)
    towards_right_lines.append(towards_right)

# ##############################################################################
# Figure setup
# ##############################################################################

linestyle_tuple = [
    ('densely dashed',        'solid'),
    ('densely dashed',        (0, (3, 2))),
    ('dashed',                (0, (7, 3))),

    ('densely dashdotted',    (0, (3, 1, 1, 1))),
    ('dashdotted',            (0, (3, 3, 1, 3))),

    ('densely dashdotdotted', (0, (5, 2, 1, 2, 1, 2))),
    ('dashdotdotted',         (0, (5, 5, 1, 5, 1, 5)))]

fig, ax = plt.subplots(1)
plt.tight_layout()
fig.set_size_inches(w=4.7747, h=3.5)

for i in range(n_heading_std):
    if i == 0:
        ax.plot(heading_angle,
                downstream_lines[i], color='b', linestyle=linestyle_tuple[i][1],
                label='downstream')
        ax.plot(heading_angle, upstream_lines[i], color='m', linestyle=linestyle_tuple[i][1],
                label='upstream')
        ax.plot(heading_angle,
                towards_left_lines[i], color='r', linestyle=linestyle_tuple[i][1],
                label='towards left')
        ax.plot(heading_angle,
                towards_right_lines[i], color='g', linestyle=linestyle_tuple[i][1],
                label='towards right')
    else:
        ax.plot(heading_angle,
                downstream_lines[i], color='b', linestyle=linestyle_tuple[i][1])
        ax.plot(heading_angle,
                upstream_lines[i], color='m', linestyle=linestyle_tuple[i][1])
        ax.plot(heading_angle,
                towards_left_lines[i], color='r', linestyle=linestyle_tuple[i][1])
        ax.plot(heading_angle,
                towards_right_lines[i], color='g', linestyle=linestyle_tuple[i][1])

ax.xaxis.set_minor_locator(MultipleLocator(np.pi / 4))
ax.xaxis.set_major_locator(MultipleLocator(base=np.pi/2))
ax.xaxis.set_major_formatter(FuncFormatter(multiple_formatter()))
ax.set_xlim([-x_max, x_max])

ax.set_ylim([0, 1.05])
# ax.yaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_major_locator(MultipleLocator(0.25))

ax.grid(True)
# ax.tick_params(axis='x')
# ax.tick_params(axis='y')

# plt.title('Semantic labels as functions of relative orientation mean and standard deviation',
#           fontdict=title_font)
plt.xlabel('Relative orientation [rad]')
plt.ylabel('Assignment confidence')

# plt.legend(loc='upper left', title='Semantic label')
plt.legend(title='Semantic label', bbox_to_anchor=(
    .5, 1.2), loc='upper center', borderaxespad=0.,  ncol=4)
#  fontsize='10', title_fontsize='10')

# ##############################################################################
# Line style legend
# ##############################################################################
std_legend_line = []
std_legend_labels = []
for i in range(n_heading_std):
    std_legend_line += ax.plot([], [], 'k', linestyle=linestyle_tuple[i][1])
    std_legend_labels.append(
        r'$\sigma_{rad}$ = ' + num_pi_formatter(sigma_heading_angle[i]))

leg = Legend(ax, std_legend_line, std_legend_labels,
             bbox_to_anchor=(1.05, .0), loc='lower left',
             title='Line style', borderaxespad=0.)
ax.add_artist(leg)


plt.savefig(
    '/home/dsp/Pictures/Matplotlib_PGFs/SemanticLabelRelativeDirection.pgf', bbox_inches='tight')
plt.show()
