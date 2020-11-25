'Interactive editor for spline evaluations'

import numpy as np
import matplotlib.pyplot as plt

from splines.AbstractSpline import Point, Points
# from splines.PolylineSpline import PolylineSpline
# from splines.HasbergSpline import HasbergSpline
from splines.CubicSpline import CubicSpline

fig, ax = plt.subplots()


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


class PathInteractor(object):
    """
    An path editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them
      'a' change aspect ratio of the figure

    """
    showverts = True
    epsilon = .5  # max distance to count as a vertex hit

    def __init__(self, figure_axes):
        self.ax = figure_axes
        canvas = self.ax.figure.canvas

        self.nodes = Points()
        self.setLemniscate()
        self._ind = None  # the active vert

        # Draw vertices and target point
        self.vertices, = self.ax.plot(self.nodes.x, self.nodes.y,
                                      marker='o',
                                      markerfacecolor='r',
                                      linestyle='None',
                                      animated=True)
        self.target_point, = ax.plot(0, -7,
                                     marker='o',
                                     markerfacecolor='b',
                                     linestyle='None',
                                     animated=True)

        # Spline container
        self.splineList = {}
        # self.addSpline(PolylineSpline('Polyline', 'k-'))
        # self.addSpline(HasbergSpline('Hasberg', 'b-'))
        self.addSpline(CubicSpline('CppCubic', 'g-'))

        self.frenetFrameList = {}
        self.constructFrenetFrames()

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event',
                           self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)

        self.canvas = canvas
        self.aspect = Aspect()
        self.relimit()

    def addSpline(self, spline):
        # Generate initial plot from nodes
        sp = spline.generatePointsFrom(self.nodes)
        self.splineList[spline], = self.ax.plot(sp.x, sp.y,
                                                spline.getFormat(),
                                                label=spline.getName(),
                                                animated=True)

    def constructFrenetFrames(self):
        target_point = self.target_point.get_xydata()
        for spline in self.splineList:
            spline.constructFrenetFrame(target_point)
            for frenet_frame in spline.frenet_frames:
                nodes = list(frenet_frame.getPoints())
                self.frenetFrameList[frenet_frame], = self.ax.plot(nodes[0],
                                                                   nodes[1])

    def updateFrenetFrames(self):
        # target_point = self.target_point.get_xydata()
        # self.frenetFrameList = {}
        # for spline in self.splineList:
        #     spline.constructFrenetFrame(target_point)
        #     for frenet_frame in spline.frenet_frames:
        #         self.frenetFrameList[frenet_frame].set_data(points)
        for k, v in self.frenetFrameList.items():
            v.remove()
        self.frenetFrameList.clear()
        self.constructFrenetFrames()

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        for spl_artist in self.splineList.values():
            self.ax.draw_artist(spl_artist)
        self.ax.draw_artist(self.vertices)
        self.ax.draw_artist(self.target_point)
        for ff_artist in self.frenetFrameList.values():
            self.ax.draw_artist(ff_artist)
        self.canvas.blit(self.ax.bbox)

    def get_ind_under_point(self, event):
        '''
        Return the index of the point closest to the event position, which is -1 for
        the target point. If no if no point is within ``self.epsilon`` to the event
        position return *None*.
        '''
        # Check target point
        d = np.sqrt((self.target_point.get_xdata()-event.xdata)**2 +
                    (self.target_point.get_ydata()-event.ydata)**2)
        if d <= self.epsilon:
            return -1

        # get the index of the vertex under nodes if within epsilon tolerance
        d = np.sqrt((self.nodes.x - event.xdata) ** 2 +
                    (self.nodes.y - event.ydata) ** 2)
        ind = d.argmin()
        if d[ind] >= self.epsilon:
            ind = None
        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        if self._ind is None:
            self.target_point.set_data(event.xdata, event.ydata)
            # Update canvas and drawings
            self.relimit()
            self.update()

        self._ind = None

    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.vertices.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        if event.key == 'a':
            self.updateAspect()
        if event.key == '+':
            self.epsilon += .1
            print("epsilon is %.1f" % self.epsilon)
        if event.key == '-':
            if self.epsilon > 0.1:
                self.epsilon -= .1
            print("epsilon is %.1f" % self.epsilon)
        if event.key == 'x':
            self.setChaos()
        if event.key == 'l':
            self.setLemniscate()
        if event.key == 'h':
            self.setHeart()

        self.update()
        self.relimit()
        self.canvas.draw()

    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        if self._ind == -1:
            self.target_point.set_data(event.xdata, event.ydata)
        else:
            x, y = event.xdata, event.ydata
            # Update point list
            self.nodes.replace(self._ind, x, y)

        # Update canvas and drawings
        self.relimit()
        self.update()

    def relimit(self):
        xmin, xmax = min(self.nodes.x), max(self.nodes.x)
        ymin, ymax = min(self.nodes.y), max(self.nodes.y)

        delta_x = (xmax - xmin) * 0.2
        delta_y = (ymax - ymin) * 0.2

        self.ax.set_xlim([xmin-delta_x, xmax+delta_x])
        self.ax.set_ylim([ymin-delta_y, ymax+delta_y])

    def updateAspect(self):
        self.aspect.switch()
        self.ax.set_aspect(self.aspect.getAspect())
        self.update()

    def update(self):
        ''' Updates the plots in the picture'''
        # Update of all structures
        # if self._ind is not None:
        # Update vertices
        self.canvas.restore_region(self.background)
        self.vertices.set_data(self.nodes.x, self.nodes.y)
        # Update splines
        for k, v in self.splineList.items():
            # Generate new points from updated nodes
            points = k.generatePointsFrom(self.nodes)
            v.set_data(points.x, points.y)
        self.updateFrenetFrames()

        self.ax.draw_artist(self.vertices)
        self.ax.draw_artist(self.target_point)
        for spl in self.splineList.keys():
            self.ax.draw_artist(self.splineList[spl])
        for ff in self.frenetFrameList.keys():
            self.ax.draw_artist(self.frenetFrameList[ff])
        self.canvas.blit(self.ax.bbox)

    def setHeart(self):
        self.nodes.x = np.array([0, 5, 4, 2, 1, 0, -1, -2, -4, -5, 0])
        self.nodes.y = np.array([-6, 0, 3, 4, 3, 0, 3, 4, 3, 0, -6])

    def setChaos(self):
        self.nodes.x = np.concatenate(
            [np.linspace(0, 20, num=10), np.linspace(10, 0, num=10)])
        self.nodes.y = np.array(
            [-4, -7, -6, -4, -6, -8, -11, -10, -5, -3, 1, 4, 5, 6, 4, 2, -2, -1, 5, 7])

    def setLemniscate(self):
        self.nodes = Points()
        alpha = 100
        t = np.linspace(0, 2*np.pi, num=30)
        self.nodes.x = np.round(alpha * np.sqrt(2) *
                                np.cos(t) / (np.sin(t)**2 + 1), 7)
        self.nodes.y = np.round(alpha * np.sqrt(2) *
                                np.cos(t) * np.sin(t) / (np.sin(t)**2 + 1), 7)


interactor = PathInteractor(ax)
ax.set_title('drag vertices to update path')
ax.set_xlabel('meter')
ax.set_ylabel('meter')
plt.autoscale(False)
plt.show()
plt.autoscale(False)
