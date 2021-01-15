import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from cycler import cycler

import corridor

from splines.AbstractSpline import Points
from splines.CubicSpline import CubicSpline


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


def confidence_ellipse(ax, mean, cov_mat, n_std=3, facecolor='none', **kwargs):
    pearson = cov_mat[0, 1]/np.sqrt(cov_mat[0, 0] * cov_mat[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)
    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov_mat[0, 0]) * n_std
    mean_x = mean[0]

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov_mat[1, 1]) * n_std
    mean_y = mean[1]

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)


def plot_cartesian_state(ax, cartesian_state, **kwargs):
    # ax.plot(cartesian_state["data"][0],
    #         cartesian_state["data"][1], 'g.', alpha=0.3)

    ax.plot(cartesian_state["mean"][0], cartesian_state["mean"][1], 'r.')
    confidence_ellipse(ax, cartesian_state["mean"], cartesian_state["cov_mat"],
                       3, facecolor='none', zorder=10, edgecolor='r', **kwargs)

    velocity_mean = cartesian_state["mean"][:2] + cartesian_state["mean"][-2:]
    velocity_cov_mat = cartesian_state["cov_mat"][-2:, -2:]
    ax.arrow(cartesian_state["mean"][0], cartesian_state["mean"][1],
             cartesian_state["mean"][2], cartesian_state["mean"][3],
             head_width=.3, head_length=.5, fc='k', hatch='o',
             length_includes_head=True)

    confidence_ellipse(ax, velocity_mean, velocity_cov_mat,
                       3, facecolor='none', zorder=10, edgecolor='b', **kwargs)


def plot_frenet_state(ax, frenet_state, **kwargs):
    # if len(frenet_state["data"]) > 0:
    #     ax.plot(frenet_state["data"][0],
    #             frenet_state["data"][1], 'g.', alpha=0.3)

    ax.plot(frenet_state["mean"][0], frenet_state["mean"][1], 'r.')
    confidence_ellipse(ax, frenet_state["mean"], frenet_state["cov_mat"],
                       3, facecolor='none', zorder=10, edgecolor='r', **kwargs)

    ax.arrow(frenet_state["mean"][0], frenet_state["mean"][1],
             frenet_state["mean"][2], frenet_state["mean"][3],
             head_width=.3, head_length=.5, fc='k', hatch='o',
             length_includes_head=True)

    velocity_mean = frenet_state["mean"][:2] + frenet_state["mean"][-2:]
    velocity_cov_mat = frenet_state["cov_mat"][-2:, -2:]
    confidence_ellipse(ax, velocity_mean, velocity_cov_mat,
                       3, facecolor='none', zorder=10, edgecolor='b', **kwargs)


def lemniscate(alpha, num_nodes):
    t_lemniscate = np.linspace(0, 2*np.pi, num=num_nodes)
    lemniscate_nodes = Points()

    lemniscate_nodes.x = alpha * np.sqrt(2) * np.cos(t_lemniscate) / \
        (np.sin(t_lemniscate)**2 + 1)

    lemniscate_nodes.y = alpha * np.sqrt(2) * np.cos(t_lemniscate) * \
        np.sin(t_lemniscate) / (np.sin(t_lemniscate)**2 + 1)

    return lemniscate_nodes


def state_sample(mean, cov_mat, num_samples):
    x, y, vx, vy = np.random.multivariate_normal(mean, cov_mat, num_samples).T

    # Calculate mean and covariance matrix of the sample data
    cart_data = np.stack((x, y, vx, vy), axis=0)
    cart_mean = np.mean(cart_data, axis=1)
    cart_cov = np.cov(cart_data)

    state = dict()
    state["data"] = cart_data
    state["mean"] = cart_mean
    state["cov_mat"] = cart_cov

    return state


def to_flat_cartesian_state(cartesian_state):
    mean = cartesian_state["mean"]
    cov_mat = cartesian_state["cov_mat"]

    state = corridor.FlatCartesianStateAndCovMat2D()
    state.x = mean[0]
    state.y = mean[1]
    state.vx = mean[2]
    state.vy = mean[3]

    state.var_x = cov_mat[0, 0]
    state.var_y = cov_mat[1, 1]
    state.var_vx = cov_mat[2, 2]
    state.var_vy = cov_mat[3, 3]
    state.cov_xy = cov_mat[0, 1]
    state.cov_xvx = cov_mat[0, 2]
    state.cov_xvy = cov_mat[0, 3]
    state.cov_yvx = cov_mat[1, 2]
    state.cov_yvy = cov_mat[1, 3]
    state.cov_vxvy = cov_mat[2, 3]

    return state


def to_frenet_state(flat_frenet_state):
    mean = np.array([flat_frenet_state.l,
                     flat_frenet_state.d,
                     flat_frenet_state.vl,
                     flat_frenet_state.vd])

    cov_mat = np.array([[flat_frenet_state.var_l, flat_frenet_state.cov_ld,
                         flat_frenet_state.cov_lvl,  flat_frenet_state.cov_lvd],
                        [flat_frenet_state.cov_ld, flat_frenet_state.var_d,
                         flat_frenet_state.cov_dvl,  flat_frenet_state.cov_dvd],
                        [flat_frenet_state.cov_lvl, flat_frenet_state.cov_dvl,
                         flat_frenet_state.var_vl,  flat_frenet_state.cov_vlvd],
                        [flat_frenet_state.cov_lvd, flat_frenet_state.cov_dvd,
                         flat_frenet_state.cov_vlvd,  flat_frenet_state.var_vd]])

    frenet_state = dict()
    frenet_state["data"] = []
    frenet_state["mean"] = mean
    frenet_state["cov_mat"] = cov_mat

    # print("UT Frenet state")
    # print(frenet_state["mean"])
    # print(frenet_state["cov_mat"])

    return frenet_state


def monte_carlo_transformation(moving_frenet_frame_assumption,
                               corridor_wrapper, cartesian_state):
    frenet_state_vector_list = list()
    # Perform transformation
    for column in cartesian_state["data"].T:
        frenet_state_vector_list.append(corridor_wrapper.to_frenet_state_vector(
            column.tolist(), moving_frenet_frame_assumption))

    # Calculate mean and cov mat
    mc_frenet_states = np.array([np.array(x)
                                 for x in frenet_state_vector_list]).T
    # mc_mean = np.mean(mc_frenet_states, axis=1)
    # mc_cov = np.cov(mc_frenet_states)

    frenet_state = dict()
    frenet_state["data"] = mc_frenet_states
    frenet_state["mean"] = np.mean(mc_frenet_states, axis=1)
    frenet_state["cov_mat"] = np.cov(mc_frenet_states)

    # print("MC Frenet state")
    # print(frenet_state["mean"])
    # print(frenet_state["cov_mat"])

    return frenet_state


def unscented_transformation(moving_frenet_frame_assumption,
                             corridor_wrapper, cartesian_state):
    flat_cartesian_state = to_flat_cartesian_state(cartesian_state)
    flat_frenet_state = corridor.ut_cartesian_frenet_transformation(
        corridor_wrapper, flat_cartesian_state, moving_frenet_frame_assumption)
    return to_frenet_state(flat_frenet_state)


def linearized_tranformation(moving_frenet_frame_assumption, corridor_wrapper, cartesian_state):
    flat_cartesian_state = to_flat_cartesian_state(cartesian_state)
    flat_frenet_state = corridor_wrapper.to_frenet_state(
        flat_cartesian_state, moving_frenet_frame_assumption)
    return to_frenet_state(flat_frenet_state)


def z_test(ground_truth, test):
    delta_mean = test["mean"] - ground_truth["mean"]
    inv_delta_cov = np.linalg.inv(
        np.linalg.inv(test["cov_mat"]) + np.linalg.inv(test["cov_mat"]))
    t1 = np.dot(delta_mean.T, inv_delta_cov)
    return np.dot(t1, delta_mean)


def generate_transformations(moving_frenet_frame_assumption,
                             corridor_wrapper, cartesian_states):
    mc_frenet_states = list()
    ut_frenet_states = list()
    ln_frenet_states = list()

    for state in cartesian_states:
        # Monte Carlo transformation
        mc_frenet_states.append(monte_carlo_transformation(
            moving_frenet_frame_assumption, corridor_wrapper, state))
        # Unscented transformation
        ut_frenet_states.append(unscented_transformation(
            moving_frenet_frame_assumption, corridor_wrapper, state))
        # Linearized transformation
        ln_frenet_states.append(linearized_tranformation(
            moving_frenet_frame_assumption, corridor_wrapper, state))

    return mc_frenet_states, ut_frenet_states, ln_frenet_states


def main():

    moving_frenet_frame_assumption = True

    # Initialize corridor from lemniscate
    lemniscate_nodes = lemniscate(30, 6)
    corridor_wrapper = corridor.CorridorWrapper(
        456, lemniscate_nodes.x.tolist(), lemniscate_nodes.y.tolist())

    corridor.TestCorridorHandle(corridor_wrapper)

    polyline_dict = corridor_wrapper.get_polylines(.1)

    # Define states
    cartesian_states = list()

    mean_1 = [12., 10., -9., -5.]
    cov_mat_1 = [[2, .5,   0,  0],
                 [.5, 1.5, 0, 0],
                 [0, 0, 0.9, 0.2],
                 [0, 0, 0.2, 1.1]]
    cartesian_states.append(state_sample(mean_1, cov_mat_1, 5000))

    mean_2 = [-28., 1., -2., 10.]
    cov_mat_2 = [[2, .5,   0,  0],
                 [.5, 1.5, 0, 0],
                 [0, 0, 1.1, 0.2],
                 [0, 0, 0.2, 1.3]]
    cartesian_states.append(state_sample(mean_2, cov_mat_2, 5000))

    mean_3 = [-24, 17, 8, -1]
    cov_mat_3 = [[2, -.5,   0,  0],
                 [-.5, 1.5, 0, 0],
                 [0, 0, 1.1, 0.2],
                 [0, 0, 0.2, 1.3]]
    cartesian_states.append(state_sample(mean_3, cov_mat_3, 5000))

    # transformations under assumption A-1
    mc_frenet_states_A1, ut_frenet_states_A1, ln_frenet_states_A1 = generate_transformations(
        False, corridor_wrapper, cartesian_states)

    # transformations under assumption A-2
    mc_frenet_states_A2, ut_frenet_states_A2, ln_frenet_states_A2 = generate_transformations(
        True, corridor_wrapper, cartesian_states)

    # State transformation evaluation
    a1_ln = []
    a1_ut = []
    a2_ln = []
    a2_ut = []
    for i in range(len(mc_frenet_states_A1)):
        # euclidean distance
        a1_ln.append(np.linalg.norm(ln_frenet_states_A1[i]['mean'] -
                                    mc_frenet_states_A1[i]['mean']))
        a1_ut.append(np.linalg.norm(ut_frenet_states_A1[i]['mean'] -
                                    mc_frenet_states_A1[i]['mean']))

        a2_ln.append(np.linalg.norm(ln_frenet_states_A2[i]['mean'] -
                                    mc_frenet_states_A2[i]['mean']))
        a2_ut.append(np.linalg.norm(ut_frenet_states_A2[i]['mean'] -
                                    mc_frenet_states_A2[i]['mean']))
        # z-test value
        a1_ln.append(z_test(mc_frenet_states_A1[i], ln_frenet_states_A1[i]))
        a1_ut.append(z_test(mc_frenet_states_A1[i], ut_frenet_states_A1[i]))

        a2_ln.append(z_test(mc_frenet_states_A2[i], ln_frenet_states_A2[i]))
        a2_ut.append(z_test(mc_frenet_states_A2[i], ut_frenet_states_A2[i]))

    table_data = [a1_ln, a1_ut, a2_ln, a2_ut]
    print(table_data)

    # PLOTTING #################################################################
    subplot_width = [0.5, 0.5]
    subplot_height = [1]
    # plt.rc('axes', prop_cycle=(cycler('color', ['m', 'g', 'b', 'c'])))
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, figure=fig)
    # , width_ratios=subplot_width,
    #  height_ratios=subplot_height)

    ax_cart = fig.add_subplot(gs[:1, :1], aspect='equal')
    ax_frenet_A1 = fig.add_subplot(gs[:1, -1], aspect='equal')
    ax_frenet_A2 = fig.add_subplot(gs[1:2, :1], aspect='equal')
    ax_table = fig.add_subplot(gs[1:2, -1])
    ax_table.axis('off')

    ax_cart.set_title('Cartesian coordinates')
    ax_frenet_A1.set_title('Frenet coordinates,\nstatic Frenet frame (A-1)')
    ax_frenet_A2.set_title('Frenet coordinates,\nmoving Frenet frame (A-2)')
    # cartesian axis labels
    ax_cart.set(xlabel='$x$ [m]')
    ax_cart.set(ylabel='$y$ [m]')

    # frenet A1 axis labels
    ax_frenet_A1.set(xlabel='$l$ [m]')
    ax_frenet_A1.set(ylabel='$d$ [m]')

    ax_cart.plot(lemniscate_nodes.x, lemniscate_nodes.y, 'k.')

    # Cartesian plot

    ax_cart.plot(polyline_dict["reference_line_x"],
                 polyline_dict["reference_line_y"], 'k-.', linewidth=1)
    ax_cart.plot(polyline_dict["left_boundary_x"],
                 polyline_dict["left_boundary_y"], color='tab:gray', linewidth=1)
    ax_cart.plot(polyline_dict["right_boundary_x"],
                 polyline_dict["right_boundary_y"], color='tab:gray', linewidth=1)

    plot_cartesian_state(ax_cart, cartesian_states[0])
    plot_cartesian_state(ax_cart, cartesian_states[1])
    plot_cartesian_state(ax_cart, cartesian_states[2])

    # Frenet plot
    ax_frenet_A1.plot([20, 130], [0, 0], 'k-.', linewidth=1)
    ax_frenet_A1.plot([20, 130], [2, 2], color='tab:gray', linewidth=1)
    ax_frenet_A1.plot([20, 130], [-2, -2], color='tab:gray', linewidth=1)

    ax_frenet_A2.plot([20, 130], [0, 0], 'k-.', linewidth=1)
    ax_frenet_A2.plot([20, 130], [2, 2], color='tab:gray', linewidth=1)
    ax_frenet_A2.plot([20, 130], [-2, -2], color='tab:gray', linewidth=1)

    # Fill A-1 plot
    for state in mc_frenet_states_A1:
        plot_frenet_state(ax_frenet_A1, state, alpha=0.5)

    for state in ut_frenet_states_A1:
        plot_frenet_state(ax_frenet_A1, state, linestyle='--', alpha=0.5)

    for state in ln_frenet_states_A1:
        plot_frenet_state(ax_frenet_A1, state, linestyle='-.', alpha=0.5)

    # Fill A-2 plot
    for state in mc_frenet_states_A2:
        plot_frenet_state(ax_frenet_A2, state, alpha=0.5)

    for state in ut_frenet_states_A2:
        plot_frenet_state(ax_frenet_A2, state, linestyle='--', alpha=0.5)

    for state in ln_frenet_states_A2:
        plot_frenet_state(ax_frenet_A2, state, linestyle='-.', alpha=0.5)

    ax_cart.set_xlim([-50, 60])
    ax_cart.set_ylim([-25, 25])

    ax_frenet_A1.set_xlim([20, 130])
    ax_frenet_A1.set_ylim([-25, 25])
    ax_frenet_A2.set_xlim([20, 130])
    ax_frenet_A2.set_ylim([-25, 25])

    # Fill table
    column_labels = ['Euclidean distance', 'Z-test values']
    row_labels = ['A-1 LT', 'A-1 UT',
                  'A-2 LT', 'A-2 UT']
    # table = ax_table.table(cellText=table_data, colLabels=column_labels,
    #                        rowLabels=row_labels, loc="center")
    # table.set_fontsize(14)

    # ax_cart.set_aspect('equal', adjustable='box')
    # ax_frenet_A1_ut.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    plt.savefig(
        '/home/dsp/Pictures/Matplotlib_PGFs/StateTransformation.pgf', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()
