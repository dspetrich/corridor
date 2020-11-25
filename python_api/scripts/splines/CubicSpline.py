import numpy as np

from .AbstractSpline import *
import corridor


class FrenetFrame:
    def __init__(self, idx, segm_arc_length, origin, tangent, normal, target_point):
        self.idx = idx
        self.segm_arc_length = segm_arc_length
        self.origin = np.array([origin['x'], origin['y']])
        self.tangent = np.array([tangent['x'], tangent['y']])
        self.normal = np.array([normal['x'], normal['y']])
        self.target_point = np.array([target_point[0][0], target_point[0][1]])

    def getPoints(self):
        return zip(self.origin, self.target_point)


class CubicSpline(AbstractSpline):
    def __init__(self, name="cubicSpline", format='g-'):
        AbstractSpline.__init__(self, name, format)
        self.frenet_frames = []

    def generatePointsFrom(self, nodes):
        self.setNodes(nodes)
        self.defineSplineParams()
        self.getPointsFromParams()
        return self.points

    def defineSplineParams(self):
        param_dict = corridor.create_spline_params(
            self.nodes.x.tolist(), self.nodes.y.tolist())

        self.parameter.h = np.diff(param_dict["arc_length"])
        self.parameter.arc_length = param_dict["arc_length"]

        self.parameter.a_x = param_dict["a_x"]
        self.parameter.b_x = param_dict["b_x"]
        self.parameter.c_x = param_dict["c_x"]
        self.parameter.d_x = param_dict["d_x"]

        self.parameter.a_y = param_dict["a_y"]
        self.parameter.b_y = param_dict["b_y"]
        self.parameter.c_y = param_dict["c_y"]
        self.parameter.d_y = param_dict["d_y"]

        self.parameter.size = len(self.parameter.h)

    def getPointsFromParams(self):
        for i in range(0, self.parameter.size):
            steps = np.linspace(0, self.parameter.h[i], num=100)
            pts = self.parameter.cartesianPosition(i, steps)
            x, y = zip(*pts)
            self.points.x = np.concatenate((self.points.x, x))
            self.points.y = np.concatenate((self.points.y, y))
        return self.getPoints()

    def constructFrenetFrame(self, target_point):
        self.frenet_frames = []
        frenet_frames = corridor.construct_frenet_frames(
            target_point[0, 0], target_point[0, 1])
        for ff in frenet_frames:
            self.frenet_frames.append(
                FrenetFrame(ff['segm_index'],
                            ff['segm_arc_length'],
                            ff['origin'],
                            ff['tangent'],
                            ff['normal'],
                            target_point))