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
        self.start_simulation(0.01, 100)

    def collision_points(self):
        """
        for i in range(len(self.plane)):
            m = min([self.plane[i].equation(x.pos) for x in self.body.particles])
            if m < 0:
                return list(filter(lambda x: abs(self.plane[i].equation(x.pos) - m) < 0.00001, self.body.particles)),i 
        return False, 1
        """
        index = 0
        for i in self.plane: 
            vert_plus = []
            vert_minus = []
            for j in self.body.particles:
                #if((j.pos[0] + self.body.centre_of_mas_pos[0])*i.unit_normal[0] + (j.pos[1] + self.body.centre_of_mas_pos[1])*i.unit_normal[1] + (j.pos[2] + self.body.centre_of_mas_pos[2])*i.unit_normal[2] -1 > 0):
                #if(j.pos[0]*i.coefficients[0] + j.pos[1]*i.coefficients[1] + j.pos[2]*i.coefficients[2] - i.coefficients[3] > 0):
                if(j.pos[0]*i.coefficients[0] + j.pos[1]*i.coefficients[1] + j.pos[2]*i.coefficients[2] - 1 > 0):
                    vert_plus.append(j)
                else:
                    vert_minus.append(j)
                if(len(vert_plus) != 0 and len(vert_minus) != 0):
                    if(len(vert_plus) <= len(vert_minus)):
                        return vert_plus, index
                    else:
                        return vert_minus, index
            index += 1
        return 0, index#"""
    
    def solve_contact(self, points, i):
        m = sum(x.mass for x in points)
        point = sum((x.mass * x.pos) for x in points) / m
        displacement = point - self.body.centre_of_mas_pos
        vel = self.body.velocity + np.cross(self.body.angular_velocity, displacement)
        relative_vel = np.dot(self.plane[i].unit_normal, vel) 
        denominator = np.cross(displacement, self.plane[i].unit_normal)
        denominator = np.matmul(denominator, np.linalg.inv(self.body.inertia_tensor).transpose())
        denominator = np.cross(denominator, displacement)
        denominator = np.dot(self.plane[i].unit_normal, denominator)
        denominator += 1 / self.body.mass
        j = -2 * relative_vel / denominator
        force = j * self.plane[i].unit_normal
        self.body.linear_momentum += force
        self.body.angular_momentum += np.cross(displacement, force)
        self.body.velocity = self.body.linear_momentum / self.body.mass
        self.body.angular_velocity = np.matmul(self.body.angular_momentum,
                                               np.linalg.inv(self.body.inertia_tensor).transpose())

    def start_simulation(self, dt, fps):
        while True:
            time.sleep(1 / fps)
            self.body.eilers_step(dt)
            c_p, number = self.collision_points()
            if c_p:
                self.solve_contact(c_p, number)
            self.body_view.update()


if __name__ == "__main__":
    cube = RigidBody(cube_generator(10, 0.2, [-10., 0., 30.]), [1., 50., 25., 100.], [0., 0., 0.], [0., 0., 2500.])
    plane = []
    plane.append(Plane([[25, 25, -8], [25, -25, 1], [-25, -25, 10], [-25, 25, 1]]))
    plane.append(Plane([[100, 40, 100], [100, 40, -100], [-100, 40,-100], [-100, 40, 100]]))
    plane.append(Plane([[100, -40, 100], [100, -40, -100], [-100, -40,-100], [-100, -40, 100]]))
    plane.append(Plane([[40, 100, 100], [40, 100, -100], [40, -100,-100], [40, 100, -100]]))
    plane.append(Plane([[-40, 100, 100], [-40, 100, -100], [-40, -100,-100], [-40, 100, -100]]))
    s = Simulator(cube, plane)
