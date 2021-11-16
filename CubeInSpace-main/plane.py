import numpy as np

#"""
class Plane:
    def __init__(self, vertices):
        self.vertices = vertices

        system_matrix = np.array([[x, y, z] for x, y, z in vertices[:3]])
        free_members = -np.ones(3)
        self.coefficients = -np.linalg.solve(system_matrix, free_members)

        self.unit_normal = self.coefficients / np.dot(self.coefficients, self.coefficients) ** 0.5

        print(self.coefficients)

    def equation(self, point):
        a, b, c = self.coefficients[0], self.coefficients[1], self.coefficients[2]
        return a * point[0] + b * point[1] + c * point[2] - 1
"""

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

class Plane:
    def __init__(self, vert):
        self.vertices = vert
        A = vert[0][1]*(vert[1][2]-vert[2][2])+vert[1][1]*(vert[2][2]-vert[0][2])+vert[2][1]*(vert[0][2]-vert[1][2])
        B = vert[0][2]*(vert[1][0]-vert[2][0])+vert[1][2]*(vert[2][0]-vert[0][0])+vert[2][2]*(vert[0][0]-vert[1][0])
        C = vert[0][0]*(vert[1][1]-vert[2][1])+vert[1][0]*(vert[2][1]-vert[0][1])+vert[2][0]*(vert[0][1]-vert[1][1])
        D = vert[0][0]*(vert[1][1]*vert[2][2]-vert[2][1]*vert[1][2])+vert[0][1]*(vert[2][1]*vert[0][2]-vert[0][1]*vert[2][2])+vert[0][2]*(vert[0][1]*vert[1][2]-vert[1][1]*vert[0][2])
        self.coefficients = [A, B, C, D]
        self.coefficients = np.asarray(self.coefficients, dtype=np.float64)
        #self.coefficients = normalize(self.coefficients)
        if D != 0: self.coefficients = self.coefficients/D
        self.unit_normal = [A, B, C]
        self.unit_normal = np.asarray(self.unit_normal, dtype=np.float64)
        self.unit_normal = normalize(self.unit_normal)
        print(self.coefficients)
        print(self.unit_normal)

    def equation(self, point):
        a, b, c, d = self.coefficients[0], self.coefficients[1], self.coefficients[2], self.coefficients[3] 
        return a * point[0] + b * point[1] + c * point[2] - 1
#"""
