"""
===========
Corridor Editor
===========

Sharing events across GUIs.

This example demonstrates a cross-GUI application using Matplotlib event
handling to interact with and modify objects on the canvas.
"""

import numpy as np
import matplotlib
from matplotlib.backend_bases import MouseButton
import matplotlib.pyplot as plt

import corridor
from base_data import Points, lemniscate


matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': '10',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.autolayout': True,
    'figure.figsize': [10, 6],
    'xtick.labelsize': 'small',
    'ytick.labelsize': 'small',
    'legend.fontsize': 'x-small',
    'legend.title_fontsize': 'small',
    'lines.linewidth': '0.5'
    # 'axes.labelsize': 'small',
})


class FrenetFrame():
    """
    Converts a frenet frame dict to a python class
    """

    def __init__(self, corridor, target_point):
        self.cartesian_point = np.array(
            [target_point[0, 0], target_point[0, 1]])
        frenet_frame_dict = corridor.get_frenet_frame_dict(
            self.cartesian_point[0], self.cartesian_point[1])
        self.origin = np.array([frenet_frame_dict["origin"]["x"],
                                frenet_frame_dict["origin"]["y"]])
        self.tangent = np.array([frenet_frame_dict["tangent"]["x"],
                                 frenet_frame_dict["tangent"]["y"]])
        self.normal = np.array([frenet_frame_dict["normal"]["x"],
                                frenet_frame_dict["normal"]["y"]])

        # Points for visualization
        self.points = Points()
        self.points.add(self.cartesian_point[0], self.cartesian_point[1])
        self.points.add(self.origin[0], self.origin[1])
        self.points.add(self.origin[0]+5*self.tangent[0],
                        self.origin[1]+5*self.tangent[1])


class Aspect():
    def __init__(self):
        self.aspectLst = ['auto', 'equal']
        self.idx = 0

    def getAspect(self):
        return self.aspectLst[self.idx]

    def switch(self):
        self.idx = (self.idx + 1) % len(self.aspectLst)
        print("New Aspect: %s" % self.getAspect())
        return self.getAspect()


class CorridorInteractor:
    """
    An path editor.

    Press 't' to toggle vertex markers on and off.  When vertex markers are on,
    they can be dragged with the mouse.
    """

    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit

    def __init__(self, fig, figure_axes):
        self.fig = fig
        self.ax = figure_axes
        self.canvas = self.ax.figure.canvas
        self._ind = None  # the active vert

        # Set a lemniscate as initial corridor nodes
        self.nodes = lemniscate(30, 17)

        self.corridor = corridor.CorridorWrapper(
            1, self.nodes.x.tolist(), self.nodes.y.tolist())
        self.plot_corridor(self.corridor)

        # Draw vertices and target point
        self.vertices, = self.ax.plot(self.nodes.x, self.nodes.y,
                                      marker='o',
                                      color='k',
                                      markerfacecolor='k',
                                      linestyle='None',
                                      animated=True,
                                      label='reference line nodes')
        self.target_point, = ax.plot(15, 5,
                                     marker='o',
                                     color='b',
                                     markerfacecolor='b',
                                     linestyle='None',
                                     animated=True,
                                     label='cartesian position')

        frenet_frame = self.constructFrenetFrame()
        self.frenet_frame_plot, = ax.plot(frenet_frame.points.x,
                                          frenet_frame.points.y,
                                          marker='.',
                                          color='m',
                                          animated=True,
                                          label='projection \& tangent')

        self.canvas.mpl_connect('draw_event', self.on_draw)

        self.canvas.mpl_connect('button_press_event',
                                self.on_button_press)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.canvas.mpl_connect('button_release_event',
                                self.on_button_release)
        self.canvas.mpl_connect('motion_notify_event',
                                self.on_mouse_move)
        self.aspect = Aspect()
        # self.relimit()

    def get_ind_under_point(self, event):
        """
        Return the index of the point closest to the event position, which is -1 for
        the target point. If no if no point is within ``self.epsilon`` to the event
        position return *None*.
        """
        # Check target point
        xyt = self.ax.transData.transform(
            (self.target_point.get_xdata(), self.target_point.get_ydata()))
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        if d <= self.epsilon:
            return -1

        if self.showverts:
            # Check nodes
            xy = np.column_stack((self.nodes.x, self.nodes.y))
            xyt = self.ax.transData.transform(xy)
            xt, yt = xyt[:, 0], xyt[:, 1]
            d = np.sqrt((xt - event.x)**2 + (yt - event.y)**2)
            ind = d.argmin()
            if d[ind] >= self.epsilon:
                return None
            return ind

    def on_draw(self, event):
        """Callback for draws."""
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.update_corridor()
        self.update_frenet_frame_plot()
        self.ax.draw_artist(self.vertices)
        self.ax.draw_artist(self.target_point)
        self.canvas.blit(self.ax.bbox)

    def on_button_press(self, event):
        """Callback for mouse button presses."""
        if (event.inaxes is None
                or event.button != MouseButton.LEFT):
            return
        self._ind = self.get_ind_under_point(event)

    def on_button_release(self, event):
        """Callback for mouse button releases."""
        if (event.button != MouseButton.LEFT):
            return
        if self._ind is None:
            self.target_point.set_data([event.xdata], [event.ydata])
            self.update()

        self._ind = None

    def on_key_press(self, event):
        """Callback for key presses."""
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.vertices.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        if event.key == 'r':
            self.nodes = lemniscate(100, 30)
            self.update()
        if event.key == 'a':
            self.updateAspect()
        self.canvas.draw()
        if event.key == 'p':
            self.fig.savefig(
                '/home/dsp/Pictures/Matplotlib_PGFs/CorridorExample.png',
                bbox_inches='tight',
                dpi=300)

    def on_mouse_move(self, event):
        """Callback for mouse movements."""
        if (self._ind is None
                or event.inaxes is None
                or event.button != MouseButton.LEFT):
            return

        if self._ind == -1:
            self.target_point.set_data([event.xdata], [event.ydata])

        elif self.showverts:
            x, y = event.xdata, event.ydata
            # Update point list
            self.nodes.replace(self._ind, x, y)

        self.update()

    def update(self):
        ''' Updates the plots in the picture'''
        # Update vertices
        self.canvas.restore_region(self.background)
        self.update_corridor()
        self.update_frenet_frame_plot()
        self.vertices.set_data(self.nodes.x, self.nodes.y)
        self.ax.draw_artist(self.vertices)
        self.ax.draw_artist(self.target_point)
        self.canvas.blit(self.ax.bbox)

    def update_corridor(self):
        self.corridor = corridor.CorridorWrapper(
            1, self.nodes.x.tolist(), self.nodes.y.tolist())
        polyline_dict = self.corridor.get_polylines(1)
        self.reference_line.set_data(polyline_dict["reference_line_x"],
                                     polyline_dict["reference_line_y"])
        self.left_boundary.set_data(polyline_dict["left_boundary_x"],
                                    polyline_dict["left_boundary_y"])
        self.right_boundary.set_data(polyline_dict["right_boundary_x"],
                                     polyline_dict["right_boundary_y"])
        self.ax.draw_artist(self.reference_line)
        self.ax.draw_artist(self.left_boundary)
        self.ax.draw_artist(self.right_boundary)

    def update_frenet_frame_plot(self):
        frenet_frame = self.constructFrenetFrame()
        self.frenet_frame_plot.set_data(frenet_frame.points.x,
                                        frenet_frame.points.y)
        self.ax.draw_artist(self.frenet_frame_plot)

    def plot_corridor(self, corridor):
        polyline_dict = corridor.get_polylines(1)
        self.reference_line, = self.ax.plot(polyline_dict["reference_line_x"],
                                            polyline_dict["reference_line_y"],
                                            'k-.',
                                            linewidth=1,
                                            animated=True,
                                            label='reference-line')
        self.left_boundary, = self.ax.plot(polyline_dict["left_boundary_x"],
                                           polyline_dict["left_boundary_y"],
                                           color='tab:red',
                                           linewidth=1,
                                           animated=True,
                                           label='left-boundary polygon')
        self.right_boundary, = self.ax.plot(polyline_dict["right_boundary_x"],
                                            polyline_dict["right_boundary_y"],
                                            color='tab:green',
                                            linewidth=1,
                                            animated=True,
                                            label='right-boundary polygon')

    def constructFrenetFrame(self):
        target_point = self.target_point.get_xydata()
        return FrenetFrame(self.corridor, target_point)

    def updateAspect(self):
        self.aspect.switch()
        self.ax.set_aspect(self.aspect.getAspect())
        # self.update()


fig, ax = plt.subplots()
interactor = CorridorInteractor(fig, ax)
# ax.set_title('drag vertices to update path')
ax.set_xlabel('meter')
ax.set_ylabel('meter')
ax.set_aspect('equal')
plt.autoscale(False)
# plt.legend(loc='upper left')

plt.legend(bbox_to_anchor=(
    .5, 1.05), loc='lower center', borderaxespad=0.,  ncol=3)


plt.show()
# fig.savefig(
#     '/home/dsp/Pictures/Matplotlib_PGFs/CorridorExample.pdf', bbox_inches='tight')
