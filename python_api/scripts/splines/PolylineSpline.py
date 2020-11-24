from .AbstractSpline import *

# a dummy spline. does nothing but giving back the points of the given nodes.
class PolylineSpline(AbstractSpline):
    
    def __init__(self, name="Polyline", format='m-x'):
        AbstractSpline.__init__(self, name, format)

    def generatePointsFrom(self, nodes, target_point):
        self.nodes = nodes
        self.points = nodes
        return self.points