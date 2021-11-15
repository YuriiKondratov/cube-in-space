import numpy as np


class Quaternion:
    def __init__(self, arr):
        self.r = arr[0]
        self.v = np.array(arr[1:])

    def __mul__(self, other):
        if type(other) == int or type(other) == float:
            return Quaternion([self.r * other, *(self.v * other)])
        r1, v1 = self.r, self.v
        r2, v2 = other.r, other.v
        res = (r1 * r2 - np.dot(v1, v2), r1 * v2 + r2 * v1 + np.cross(v1, v2))
        return Quaternion([res[0], *res[1]])

    def __add__(self, other):
        return Quaternion([self.r + other.r, *(self.v + other.v)])

    def __str__(self):
        return str([self.r] + [*self.v])

    def normalize(self):
        n = (self.r ** 2 + np.dot(self.v, self.v)) ** 0.5
        return Quaternion([self.r / n, *(self.v / n)])

    def qtom(self):
        q = [self.r, *self.v]
        return np.array([[2 * (q[0] * q[0] + q[1] * q[1]) - 1,
                          2 * (q[1] * q[2] - q[0] * q[3]),
                          2 * (q[1] * q[3] + q[0] * q[2])],
                         [2 * (q[1] * q[2] + q[0] * q[3]),
                          2 * (q[0] * q[0] + q[2] * q[2]) - 1,
                          2 * (q[2] * q[3] - q[0] * q[1])],
                         [2 * (q[1] * q[3] - q[0] * q[2]),
                          2 * (q[2] * q[3] + q[0] * q[1]),
                          2 * (q[0] * q[0] + q[3] * q[3]) - 1]])
