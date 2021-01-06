import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cycler import cycler

import corridor

from splines.AbstractSpline import Points
from splines.CubicSpline import CubicSpline


def lemniscate(alpha, num_nodes):
    t_lemniscate = np.linspace(0, 2*np.pi, num=num_nodes)
    lemniscate_nodes = Points()

    lemniscate_nodes.x = alpha * np.sqrt(2) * np.cos(t_lemniscate) / \
        (np.sin(t_lemniscate)**2 + 1)

    lemniscate_nodes.y = alpha * np.sqrt(2) * np.cos(t_lemniscate) * \
        np.sin(t_lemniscate) / (np.sin(t_lemniscate)**2 + 1)

    return lemniscate_nodes


def plot(lemniscate_nodes):
    fig, ax = plt.subplots()
    ax.plot(lemniscate_nodes.x, lemniscate_nodes.y)

    # spline object for visualization purposes
    spline = CubicSpline('CppCubic', 'g-')
    natural_sp = spline.generatePointsFrom(lemniscate_nodes)
    ax.plot(natural_sp.x, natural_sp.y)

    plt.show()


def main():
    lemniscate_nodes = lemniscate(100, 6)
    corridor_wrapper = corridor.CorridorWrapper(
        456, lemniscate_nodes.x.tolist(), lemniscate_nodes.y.tolist())

    corridor.TestCorridorHandle(corridor_wrapper)
    plot(lemniscate_nodes)


if __name__ == "__main__":
    main()
