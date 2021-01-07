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


def plot_cartesian_state(ax, cartesian_state):
    # ax.plot(cartesian_state["data"][0],
    #         cartesian_state["data"][1], 'g.', alpha=0.3)

    ax.plot(cartesian_state["mean"][0], cartesian_state["mean"][1], 'ro')
    confidence_ellipse(ax, cartesian_state["mean"], cartesian_state["cov_mat"],
                       3, facecolor='none', zorder=10, edgecolor='r')

    ax.arrow(cartesian_state["mean"][0], cartesian_state["mean"][1],
             cartesian_state["mean"][2], cartesian_state["mean"][3],
             head_width=.3, head_length=.5, fc='b', hatch='o')

    velocity_mean = cartesian_state["mean"][:2] + cartesian_state["mean"][-2:]
    velocity_cov_mat = cartesian_state["cov_mat"][-2:, -2:]
    confidence_ellipse(ax, velocity_mean, velocity_cov_mat,
                       3, facecolor='none', zorder=10, edgecolor='b')


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


def flat_cartesian_state(cartesian_state):
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


def main():
    lemniscate_nodes = lemniscate(30, 6)
    corridor_wrapper = corridor.CorridorWrapper(
        456, lemniscate_nodes.x.tolist(), lemniscate_nodes.y.tolist())

    corridor.TestCorridorHandle(corridor_wrapper)

    polyline_dict = corridor_wrapper.get_polylines(.1)

    # Define sample states
    mean_1 = [-24., 17., 7., -3.]
    cov_mat_1 = [[2, -.5,   0,  0],
                 [-.5, 1.5, 0, 0],
                 [0, 0, 1.1, 0.2],
                 [0, 0, 0.2, 1.3]]
    cartesian_state_1 = state_sample(mean_1, cov_mat_1, 5000)

    mean_2 = [-28., 0., 0., 10.]
    cov_mat_2 = [[2, .5,   0,  0],
                 [.5, 1.5, 0, 0],
                 [0, 0, 1.1, 0.2],
                 [0, 0, 0.2, 1.3]]
    cartesian_state_2 = state_sample(mean_2, cov_mat_2, 5000)

    mean_3 = [17., 8., -10., 0.]
    cov_mat_3 = [[2, .5,   0,  0],
                 [.5, 1.5, 0, 0],
                 [0, 0, 1.1, 0.2],
                 [0, 0, 0.2, 1.3]]
    cartesian_state_3 = state_sample(mean_3, cov_mat_3, 5000)

    flat_cartesian_state_1 = flat_cartesian_state(cartesian_state_1)
    flat_cartesian_state_2 = flat_cartesian_state(cartesian_state_2)
    flat_cartesian_state_3 = flat_cartesian_state(cartesian_state_3)

    # Monte Carlo transformation

    # State 1
    frenet_state_vector_list = list()
    for column in cartesian_state_1["data"].T:
        frenet_state_vector_list.append(corridor_wrapper.to_frenet_state_vector(
            column.tolist()))

    mc_frenet_states = np.array([np.array(x)
                                 for x in frenet_state_vector_list]).T

    mc_mean = np.mean(mc_frenet_states, axis=1)
    mc_cov = np.cov(mc_frenet_states)
    print(mc_mean)
    print(mc_cov)

    print('before flat_cartesian_state_1')
    print(flat_cartesian_state_1.vx)
    print(flat_cartesian_state_1.vy)
    # Unscented transformation
    ut_frenet_state = corridor.ut_cartesian_frenet_transformation(
        corridor_wrapper, flat_cartesian_state_1)

    print('after flat_cartesian_state_1')
    print(flat_cartesian_state_1.vx)
    print(flat_cartesian_state_1.vy)

    print(ut_frenet_state.l)
    print(ut_frenet_state.d)
    print(ut_frenet_state.vl)
    print(ut_frenet_state.vd)

    print(ut_frenet_state.var_l)
    print(ut_frenet_state.var_d)

    # print(frenet_state_vector_list)

    # frenet_state_list = corridor_wrapper.to_frenet_state_vector()

    # PLOTTING ########################################
    fig, ax = plt.subplots()
    # ax.plot(lemniscate_nodes.x, lemniscate_nodes.y)

    # spline object for visualization purposes
    # spline = CubicSpline('CppCubic', 'g-')
    # natural_sp = spline.generatePointsFrom(lemniscate_nodes)
    # ax.plot(natural_sp.x, natural_sp.y)

    ax.plot(polyline_dict["reference_line_x"],
            polyline_dict["reference_line_y"], 'k-.')
    ax.plot(polyline_dict["left_boundary_x"],
            polyline_dict["left_boundary_y"], color='tab:gray')
    ax.plot(polyline_dict["right_boundary_x"],
            polyline_dict["right_boundary_y"], color='tab:gray')

    plot_cartesian_state(ax, cartesian_state_1)
    plot_cartesian_state(ax, cartesian_state_2)
    plot_cartesian_state(ax, cartesian_state_3)

    plt.show()


if __name__ == "__main__":
    main()
