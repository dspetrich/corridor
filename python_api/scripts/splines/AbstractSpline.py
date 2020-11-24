from ast import NodeVisitor
import numpy as np

class Point():
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Points():
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])

    def __add__(self, other):
        self.x = np.concatenate((self.x, other.x))
        self.y = np.concatenate((self.y, other.y))
        return self

    def vertices(self):
        return zip(self.x, self.y)

    def __str__(self):
        return ('x: ' + str(self.x) + '\ny: ' + str(self.y))
    
    def add(self, x, y):
        self.x = np.append(self.x, x) 
        self.y = np.append(self.y, y)
        
    def replace(self, idx, x, y):
        self.x[idx] = x
        self.y[idx] = y
        

class SplineParameter():
    def __init__(self):
        self.h = []
        self.arc_length = []
        self.size = 0
        self.a_x = []
        self.a_y = []
        self.b_x = []
        self.b_y = [] 
        self.c_x = [] 
        self.c_y = []
        self.d_x = []
        self.d_y = []

    def __add__(self, other):
        self.h += other.h
        self.arc_length += other.arc_length
        self.size += other.size
        self.a_x += other.a_x
        self.a_y += other.a_y
        self.b_x += other.b_x
        self.b_y += other.b_y
        self.c_x += other.c_x
        self.c_y += other.c_y
        self.d_x += other.d_x
        self.d_y += other.d_y
        return self

    def printParams(self):
        print(self.size)
        print(self.arc_length)
        print(self.h)
        print((self.a_x, self.a_y))
        print((self.b_x, self.b_y))
        print((self.c_x, self.c_y))
        print((self.d_x, self.d_y))
        
    def getArclengthAtIdx(self, idx):
        return self.arc_length[idx]
        
    def cartesianPosition(self, idx, local_s):
        x_coef = [self.d_x[idx], self.c_x[idx],
                  self.b_x[idx], self.a_x[idx]]
        y_coef = [self.d_y[idx], self.c_y[idx],
                  self.b_y[idx], self.a_y[idx]]
        x = np.polyval(x_coef, local_s)
        y = np.polyval(y_coef, local_s)
        # print("evaluatePosition:", x, y)
        return np.column_stack((x, y))
    
    def cartesianTangent(self, idx, local_s):
        x_coef = [3*self.d_x[idx], 2*self.c_x[idx], self.b_x[idx]]
        y_coef = [3*self.d_y[idx], 2*self.c_y[idx], self.b_y[idx]]
        x = np.polyval(x_coef, local_s)
        y = np.polyval(y_coef, local_s)
        print("cartesianTangent: ", x, y)
        return np.column_stack((x, y))
    
    def cartesianNormal(self, idx, local_s):
        x_coef = [-3*self.d_y[idx], -2*self.c_y[idx], -self.b_y[idx]]
        y_coef = [3*self.d_x[idx], 2*self.c_x[idx], self.b_x[idx]]
        x = np.polyval(x_coef, local_s)
        y = np.polyval(y_coef, local_s)
        # print("evaluateNormal:", x, y)
        return np.column_stack((x, y))

    def curvatureVector(self, idx, local_s):
        x_coef = [6*self.d_x[idx], 2*self.c_x[idx]]
        y_coef = [6*self.d_y[idx], 2*self.c_y[idx]]
        x = np.polyval(x_coef, local_s)
        y = np.polyval(y_coef, local_s)
        # print("evaluateCurvature: {0}, {1}", x, y)
        return np.column_stack((x, y))

# spline interface
# a spline displayed by spline_editor.py must inherit this class. (or at least implement the same methods)
class AbstractSpline():
    # define global variables
    def __init__(self, name, formating='k-'):
        self.type = 'spline'
        self.name = name
        self.format = formating
        self.nodes = Points()
        self.points = Points()
        self.parameter = SplineParameter()
        self.visible = True

    def getName(self):
        return self.name

    def getNodes(self):
        return self.nodes

    def setNodes(self, nodes):
        self.clear()
        self.nodes = nodes

    def getPoints(self):
        return self.points

    def getFormat(self):
        return self.format

    def setVisible(self, visible):
        self.visible = visible

    def isVisible(self):
        return self.visible

    def generatePointsFrom(self, x, y, target_point):
        return self.getPoints()

    def clear(self):
        self.nodes = Points()
        self.points = Points()
        
    def findClosestNodeIdx(self, point):
        x, y = point[0]
        d = np.sqrt((self.nodes.x - x) ** 2 + 
                    (self.nodes.y - y) ** 2)
        return d.argmin()