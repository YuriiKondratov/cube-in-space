import numpy as np

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
