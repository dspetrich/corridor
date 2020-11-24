import numpy as np
import matplotlib.pyplot as plt

import corridor

from splines.AbstractSpline import Points
from splines.CubicSpline import CubicSpline


def main():
    alpha = 10

    t_lemniscate = np.linspace(0, 2*np.pi, num=1000)
    lemniscate_nodes = Points()
    lemniscate_nodes.x = alpha * np.sqrt(2) * np.cos(t_lemniscate) / \
        (np.sin(t_lemniscate)**2 + 1)
    lemniscate_nodes.y = alpha * np.sqrt(2) * np.cos(t_lemniscate) * \
        np.sin(t_lemniscate) / (np.sin(t_lemniscate)**2 + 1)

    # Knots of the Lemniscate
    plt.plot(lemniscate_nodes.x, lemniscate_nodes.y, 'r-')

    colormap = ('m', 'c', 'b', 'c')
    n_nodes = (7, 11, 21, 33)

    for k, item in enumerate(n_nodes):
        t_nodes = np.linspace(0, 2*np.pi, num=item)
        nodes = Points()
        nodes.x = alpha * np.sqrt(2) * np.cos(t_nodes) / \
            (np.sin(t_nodes)**2 + 1)
        nodes.y = alpha * np.sqrt(2) * np.cos(t_nodes) * \
            np.sin(t_nodes) / (np.sin(t_nodes)**2 + 1)

        # Knots of the Lemniscate
        plt.plot(nodes.x, nodes.y, 'x--', color=colormap[k], linewidth=0.5)

        # Cubic spline
        spline = CubicSpline('CppCubic', 'g-')
        sp = spline.generatePointsFrom(nodes)

        plt.plot(sp.x, sp.y, '-', color=colormap[k])

    plt.show()


if __name__ == "__main__":
    main()
