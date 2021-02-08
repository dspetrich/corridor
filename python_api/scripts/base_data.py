import numpy as np


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

    @classmethod
    def set(cls, x, y):
        cls.x = np.array([x])
        cls.y = np.array([y])
        return cls


def lemniscate(alpha, num_nodes):
    t_lemniscate = np.linspace(0, 2*np.pi, num=num_nodes)
    lemniscate_nodes = Points()

    lemniscate_nodes.x = alpha * np.sqrt(2) * np.cos(t_lemniscate) / \
        (np.sin(t_lemniscate)**2 + 1)

    lemniscate_nodes.y = alpha * np.sqrt(2) * np.cos(t_lemniscate) * \
        np.sin(t_lemniscate) / (np.sin(t_lemniscate)**2 + 1)

    return lemniscate_nodes
