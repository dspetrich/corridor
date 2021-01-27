import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cycler import cycler

import corridor
from base_data import Points
# from splines.AbstractSpline import Points
# from splines.CubicSpline import CubicSpline

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


def main():
    # Figure
    subplot_width = [3, 3]
    subplot_height = [2, 2]
    plt.rc('axes', prop_cycle=(cycler('color', ['m', 'g', 'b', 'c'])))
    fig = plt.figure()
    gs = gridspec.GridSpec(
        2, 2, width_ratios=subplot_width, height_ratios=subplot_height)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    # lemniscate definition
    alpha = 10

    t_lemniscate = np.linspace(0, 2*np.pi, num=1000)
    lemniscate_nodes = Points()

    lemniscate_nodes.x = alpha * np.sqrt(2) * np.cos(t_lemniscate) / \
        (np.sin(t_lemniscate)**2 + 1)

    lemniscate_nodes.y = alpha * np.sqrt(2) * np.cos(t_lemniscate) * \
        np.sin(t_lemniscate) / (np.sin(t_lemniscate)**2 + 1)

    # Plot lemniscate
    ax1.plot(lemniscate_nodes.x, lemniscate_nodes.y, 'r-')  # ,  linewidth=0.5)
    ax2.plot(lemniscate_nodes.x, lemniscate_nodes.y, 'r-')  # ,  linewidth=0.5)
    ax3.plot(lemniscate_nodes.x, lemniscate_nodes.y, 'r-')  # ,  linewidth=0.5)

    length_lemniscate = 0
    for i in range(len(t_lemniscate)-1):
        delta_x = lemniscate_nodes.x[i+1] - lemniscate_nodes.x[i]
        delta_y = lemniscate_nodes.y[i+1] - lemniscate_nodes.y[i]
        length_lemniscate += np.sqrt(delta_x**2 + delta_y**2)

    # Tangent vector to the first and last knot
    first_tangent = [0, 1]
    last_tangent = [0, 1]

    # Draw tangent at first and last point
    ax1.plot(lemniscate_nodes.x[0], lemniscate_nodes.y[0], 'k.')
    ax1.arrow(lemniscate_nodes.x[0], lemniscate_nodes.y[0],
              0, 3, head_width=.3, head_length=.5, fc='k', hatch='o')
    ax3.plot(lemniscate_nodes.x[0], lemniscate_nodes.y[0], 'k.')
    ax3.arrow(lemniscate_nodes.x[0], lemniscate_nodes.y[0],
              0, 3, head_width=.3, head_length=.5, fc='k', hatch='o')

    plotable_approximation = [7, 13, 23]

    # Storage for arc-length approximation error
    natural_arc_length_list = []
    clamped_arc_length_list = []

    natural_spline_line_handles = []
    natural_spline_labels = []
    clamped_spline_line_handles = []
    clamped_spline_labels = []

    # Range of node samples for curve approximation
    n_nodes = range(4, 40, 1)
    for k, item in enumerate(n_nodes):
        t_nodes = np.linspace(0, 2*np.pi, num=item)
        nodes = Points()
        nodes.x = alpha * np.sqrt(2) * np.cos(t_nodes) / \
            (np.sin(t_nodes)**2 + 1)
        nodes.y = alpha * np.sqrt(2) * np.cos(t_nodes) * \
            np.sin(t_nodes) / (np.sin(t_nodes)**2 + 1)

        # Natural cubic spline
        natural_sp = corridor.CubicSplineWrapper(
            1, nodes.x.tolist(), nodes.y.tolist())
        natural_arc_length_list.append(natural_sp.total_length())

        # Clamped cubic spline
        clamped_sp = corridor.CubicSplineWrapper(
            2, nodes.x.tolist(), nodes.y.tolist(), first_tangent, last_tangent)
        clamped_arc_length_list.append(clamped_sp.total_length())

        if item in plotable_approximation:
            # Natural spline
            node_plot = ax2.plot(nodes.x, nodes.y, '.-.', linewidth=0.3)
            polyline_dict = natural_sp.get_polyline(0.5)
            natural_spline_line_handles += ax2.plot(polyline_dict['x'], polyline_dict['y'], '-',
                                                    color=node_plot[0].get_color(), linewidth=.7)
            natural_spline_labels.append('n = %i' % item)
            # Clamped spline
            ax3.plot(nodes.x, nodes.y, '.-.',
                     color=node_plot[0].get_color(), linewidth=0.3)
            polyline_dict = clamped_sp.get_polyline(0.5)
            clamped_spline_line_handles += ax3.plot(polyline_dict['x'], polyline_dict['y'], '-',
                                                    color=node_plot[0].get_color(), linewidth=0.7)
            clamped_spline_labels.append('n = %i' % item)

    # Plot arc length approximation error for both cubic spline types
    natural_arc_length_list -= length_lemniscate
    natural_arc_length_list /= length_lemniscate
    clamped_arc_length_list -= length_lemniscate
    clamped_arc_length_list /= length_lemniscate
    line_handles = []
    line_handles += ax4.plot(n_nodes, abs(natural_arc_length_list), '.--')
    line_handles += ax4.plot(n_nodes, abs(clamped_arc_length_list), '.--')
    ax4.legend(line_handles, ['Natural spline',
                              'Clamped spline'], loc='upper right')

    # #####################################################################
    # Figure setup
    # #####################################################################
    ax1.set_title('Lemniscate of Bernoulli')
    ax2.set_title('Natural cubic splines')
    ax3.set_title('Clamped cubic splines')
    ax4.set_title(
        'Arc-length difference between\nlemniscate and cubic splines')

    # asp = np.diff(ax2.get_xlim())[0]/np.diff(ax2.get_ylim())[0]
    # ax4.set_aspect(asp)

    # ax1.spines["left"].set_position('center')
    # ax1.spines["bottom"].set_position('center')

    ax1.axis('equal')
    ax2.axis('equal')
    ax3.axis('equal')

    # axis labels
    ax1.set(xlabel='x [m]')
    ax1.set(ylabel='y [m]')
    ax2.set(xlabel='x [m]')
    ax2.set(ylabel='y [m]')
    ax3.set(xlabel='x [m]')
    ax3.set(ylabel='y [m]')
    ax4.set(xlabel='Number of support points')
    ax4.set(ylabel='Relative Difference [\%]')

    # legends
    ax2.legend(natural_spline_line_handles,
               natural_spline_labels, loc='upper left', ncol=4)
    ax3.legend(clamped_spline_line_handles,
               clamped_spline_labels, loc='upper left', ncol=4)

    plt.savefig(
        '/home/dsp/Pictures/Matplotlib_PGFs/LemniscateApproximation.pgf', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()
