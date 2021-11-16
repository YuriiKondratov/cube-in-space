import vpython
import numpy as np
from quaterion import Quaternion
from copy import deepcopy

g = 9.81


def star_matrix(vc):
    return np.array([[0, -vc[2], vc[1]],
                     [vc[2], 0, -vc[0]],
                     [-vc[1], vc[0], 0]])


class Particle:
    def __init__(self, p, m):
        self.pos = np.array(p)  # vector
        self.mass = m  # scalar


class RigidBody:
    def __init__(self, ps, r_q, l_m, a_m):
        vpython.scene.width = 1520
        vpython.scene.height = 725
        
        self.particles = ps
        self.mass = sum(x.mass for x in self.particles)

        # state
        self.centre_of_mas_pos = sum((x.pos * x.mass) for x in self.particles) / self.mass
        self.rotation_quaternion = Quaternion(r_q)
        self.linear_momentum = np.array(l_m)  # vector
        self.angular_momentum = np.array(a_m)  # vector

        self.displacements = deepcopy(ps)
        for x in self.displacements:
            x.pos = x.pos - self.centre_of_mas_pos

        self.body_inertia_tensor = sum(x.mass * (np.inner(x.pos, x.pos) * np.eye(3) - np.outer(x.pos, x.pos))
                                       for x in self.displacements)

        # auxiliary
        self.rotation_matrix = self.rotation_quaternion.normalize().qtom()
        self.velocity = self.linear_momentum / self.mass
        self.inertia_tensor = np.matmul(np.matmul(self.rotation_matrix, self.body_inertia_tensor),
                                        self.rotation_matrix.transpose())  # 3x3 matrix
        self.angular_velocity = np.matmul(self.angular_momentum, np.linalg.inv(self.inertia_tensor).transpose())  # vector
        self.total_force = self.mass * np.array([0., 0., -g])
        self.torque = 0

    def derivative(self):
        self.velocity = self.linear_momentum / self.mass
        self.inertia_tensor = np.matmul(np.matmul(self.rotation_matrix, self.body_inertia_tensor),
                                        self.rotation_matrix.transpose())  # 3x3 matrix
        self.angular_velocity = np.matmul(self.angular_momentum, np.linalg.inv(self.inertia_tensor).transpose())  # vector

    def eilers_step(self, dt):
        # full body state update
        self.centre_of_mas_pos += self.linear_momentum / self.mass * dt
        self.rotation_quaternion += Quaternion([0, *self.angular_velocity]) * self.rotation_quaternion * 0.5 * dt
        self.linear_momentum += self.total_force * dt
        self.angular_momentum += self.torque * dt

        # particles positions update
        self.rotation_matrix = self.rotation_quaternion.normalize().qtom()
        for particle, displacement in zip(self.particles, self.displacements):
            particle.pos = self.centre_of_mas_pos + np.matmul(displacement.pos, self.rotation_matrix.transpose())
        self.derivative()
