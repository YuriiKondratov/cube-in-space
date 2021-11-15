import numpy as np
from body import RigidBody, Particle
from view import RigidBodyView, PlaneView
from plane import Plane
import time


def cube_generator(length, r, c):
    magica = 0
    if length % 2 == 0:
        magica = r
    cube = [Particle(np.array([c[0] + r * 2 * (length / 2) - r, c[1] + r * 2 * (length / 2) - r,
                               c[2] - r * 2 * (length // 2) + r * 2 * i + magica]), 1) for i in range(length)]
    cube.extend([Particle(np.array([c[0] + r * 2 * (length / 2) - r, c[1] - r * 2 * (length / 2) + r,
                                    c[2] - r * 2 * (length // 2) + r * 2 * i + magica]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length / 2) + r, c[1] - r * 2 * (length / 2) + r,
                                    c[2] - r * 2 * (length // 2) + r * 2 * i + magica]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length / 2) + r, c[1] + r * 2 * (length / 2) - r,
                                    c[2] - r * 2 * (length // 2) + r * 2 * i + magica]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length // 2) + r * 2 * i + magica, c[1] + r * 2 * (length / 2) - r,
                                    c[2] + r * 2 * (length / 2) - r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length // 2) + r * 2 * i + magica, c[1] - r * 2 * (length / 2) + r,
                                    c[2] + r * 2 * (length / 2) - r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length // 2) + r * 2 * i + magica, c[1] - r * 2 * (length / 2) + r,
                                    c[2] - r * 2 * (length / 2) + r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length // 2) + r * 2 * i + magica, c[1] + r * 2 * (length / 2) - r,
                                    c[2] - r * 2 * (length / 2) + r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] + r * 2 * (length / 2) - r, c[1] - r * 2 * (length // 2) + r * 2 * i + magica,
                                    c[2] - r * 2 * (length / 2) + r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] + r * 2 * (length / 2) - r, c[1] - r * 2 * (length // 2) + r * 2 * i + magica,
                                    c[2] + r * 2 * (length / 2) - r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length / 2) + r, c[1] - r * 2 * (length // 2) + r * 2 * i + magica,
                                    c[2] + r * 2 * (length / 2) - r]), 1) for i in range(length)])
    cube.extend([Particle(np.array([c[0] - r * 2 * (length / 2) + r, c[1] - r * 2 * (length // 2) + r * 2 * i + magica,
                                    c[2] - r * 2 * (length / 2) + r]), 1) for i in range(length)])
    return cube


class Simulator:
    def __init__(self, rb: RigidBody, p: Plane):
        self.body = rb
        self.plane = p
        self.body_view: RigidBodyView = RigidBodyView(rb)
        self.plane_view: PlaneView = PlaneView(p)
        self.start_simulation(0.01)

    def collision_points(self):
        m = min([self.plane.equation(x.pos) for x in self.body.particles])
        if m < 0:
            return list(filter(lambda x: abs(self.plane.equation(x.pos) - m) < 0.00001, self.body.particles))
        return False

    def solve_contact(self, points):
        m = sum(x.mass for x in points)
        point = sum((x.mass * x.pos) for x in points) / m
        displacement = point - self.body.centre_of_mas_pos
        vel = self.body.velocity + np.cross(self.body.angular_velocity, displacement)
        relative_vel = np.dot(self.plane.unit_normal, vel)
        denominator = np.cross(displacement, self.plane.unit_normal)
        denominator = np.matmul(denominator, np.linalg.inv(self.body.inertia_tensor).transpose())
        denominator = np.cross(denominator, displacement)
        denominator = np.dot(self.plane.unit_normal, denominator)
        denominator += 1 / self.body.mass
        j = -2 * relative_vel / denominator
        force = j * self.plane.unit_normal
        self.body.linear_momentum += force
        self.body.angular_momentum += np.cross(displacement, force)
        self.body.velocity = self.body.linear_momentum / self.body.mass
        self.body.angular_velocity = np.matmul(self.body.angular_momentum,
                                               np.linalg.inv(self.body.inertia_tensor).transpose())

    def start_simulation(self, dt):
        while True:
            time.sleep(dt)
            self.body.eilers_step(dt)
            c_p = self.collision_points()
            if c_p:
                self.solve_contact(c_p)
            self.body_view.update()


if __name__ == "__main__":
    cube = RigidBody(cube_generator(10, 0.2, [0., 0., 15.]), [1., 50., 25., 100.], [0., 0., 0.], [0., 0., 2500.])
    plane = Plane([[25, 25, -8], [25, -25, 1], [-25, -25, 10], [-25, 25, 1]])
    s = Simulator(cube, plane)
