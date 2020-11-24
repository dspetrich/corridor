# Cubic spline, based on Paper von Hasberg 

import numpy as np
from math import *

from AbstractSpline import *

# spline implementiert nach Hasberg, siehe Doktorarbeit von Hasberg
class HasbergSpline(AbstractSpline):
    
    def __init__(self, name="hasberg", format='r-'):
        AbstractSpline.__init__(self, name, format)

    def generatePointsFrom(self, nodes, target_point):
        self.setNodes(nodes)
        params = self.getParams()
        return self.getPointsFromParams(params)
    
    def getPointsFromParams(self, params):
        for i in range(0, params.size):
            steps = np.linspace(0, params.h[i], num=100)
           
            x_coef = [params.d_x[i], params.c_x[i], params.b_x[i], params.a_x[i]]
            y_coef = [params.d_y[i], params.c_y[i], params.b_y[i], params.a_y[i]]
            
            new_x_points = np.polyval(x_coef, steps)#[0]
            new_y_points = np.polyval(y_coef, steps)#[0]
            
            self.points.x = np.concatenate((self.points.x, new_x_points))
            self.points.y = np.concatenate((self.points.y, new_y_points))
            
        return self.getPoints()
    
    def getParams(self):        
        h = self.get_h()
        M = self.get_M(h)
        L = self.get_L(h)
        M_inv = self.inverse(M)
        m = self.get_m(M_inv, L)        
        params = self.get_abcd(m[0], m[1], h)        
        return params
        
        
    def get_h(self):
        h = []
        x_nodes = self.getNodes().x
        y_nodes = self.getNodes().y
        
        if len(x_nodes) is not len(y_nodes):
            print("HasbergSpline: get_h(): x and y have not the same length.")
            return
        
        for i in range(0, len(x_nodes) - 1):
            dx = x_nodes[i+1] - x_nodes[i]
            dy = y_nodes[i+1] - y_nodes[i]
            h.append(sqrt(dx**2 + dy**2))
            
        return h
            
    def get_M(self, h):
        m_max = len(h) + 1        
        M = np.zeros((m_max, m_max))
        
        for i in range(1, m_max-1):
            M[i][i-1] = h[i-1]       
        for i in range(1, m_max-1):
            M[i][i+1] = h[i]       
        for i in range(1, m_max-1):
            M[i][i] = 2*(h[i-1] + h[i])
    
        M[0][0] = 1
        # M[0][1] = -1
        # M[m_max-1][m_max-2] = 1
        M[m_max-1][m_max-1] = -1

        M[0][0] = 1.
        M[m_max-1][m_max-1] = 1.

        return M
            
    def get_L(self, h):        
        m_max = len(h) + 1
        L = np.zeros((m_max, m_max))
        
        for i in range(1, m_max-1):
            L[i][i-1] = 6. / h[i-1]
            
        for i in range(1, m_max-1):
            L[i][i+1] = 6. / h[i]
                         
        for i in range(1, m_max-1):
            L[i][i] = (-6./h[i-1] - 6./h[i])      
        
        return L 
            
            
    def get_m(self, M_inv, L):
        x_nodes = self.getNodes().x
        y_nodes = self.getNodes().y
        tmp = np.dot(M_inv, L)

        m_x = np.dot(tmp, x_nodes)
        m_y = np.dot(tmp, y_nodes)

        return (m_x, m_y)
    
    
    def get_abcd(self, m_x, m_y, h):
        params = SplineParameter()
        
        params.size = len(h) 
        params.h = h
        
        x_nodes = self.getNodes().x
        y_nodes = self.getNodes().y
        
        for i in range(0, params.size):
            params.a_x.append(x_nodes[i])
            params.a_y.append(y_nodes[i])

            params.b_x.append( (x_nodes[i+1] - x_nodes[i]) / params.h[i] - 1./6. * params.h[i] * (m_x[i+1] + 2 * m_x[i]))
            params.b_y.append( (y_nodes[i+1] - y_nodes[i]) / params.h[i] - 1./6. * params.h[i] * (m_y[i+1] + 2 * m_y[i]))
        
            params.c_x.append(m_x[i] / 2.)
            params.c_y.append(m_y[i] / 2.)
            
            params.d_x.append( (m_x[i+1] - m_x[i]) / (6 * params.h[i]) )
            params.d_y.append( (m_y[i+1] - m_y[i]) / (6 * params.h[i]) )

        return params
     
            
    def inverse(self, A):
        return np.linalg.inv(A)   
    
    #winkel zwischen den ableitungen
    def getDerivativeAngles(self, nodes):
        self.setNodes((nodes.x, nodes.y))
        params = self.getParams()
        
        psis = []
        x_coef = []
        y_coef = []
        
        for i in range(0, params.size):

            x_coef = [3*params.d_x[i], 2*params.c_x[i], params.b_x[i]]
            y_coef = [3*params.d_y[i], 2*params.c_y[i], params.b_y[i]]
                                
            yval = np.polyval(y_coef, 0.)
            xval = np.polyval(x_coef, 0.)
            psis.append(atan2(yval, xval))
        
        yval = np.polyval(y_coef, params.h[-1])
        xval = np.polyval(x_coef, params.h[-1])
        psis.append(atan2(yval, xval))  
        
        return psis
            
