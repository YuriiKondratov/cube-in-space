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
        ind = []
        vert_plus = []
        vert_minus = []
        index = 0
        for i in self.plane: 
            vert_plus.append([])
            vert_minus.append([])
            for j in self.body.particles:
                if(j.pos[0]*i.coefficients[0] + j.pos[1]*i.coefficients[1] + j.pos[2]*i.coefficients[2] - 1 > 0):
                    vert_plus[index].append(j)
                else:
                    vert_minus[index].append(j)
            if(len(vert_plus[index]) != 0 and len(vert_minus[index]) != 0):
                ind.append(index)
            index += 1
        res = []
        if(len(ind) != 0):
            for i in range(len(vert_plus)):
                if vert_plus[i] != [] and vert_minus[i] != []:
                    if len(vert_plus[i]) <= len(vert_minus[i]):
                        res.append(vert_plus[i])
                    else:
                        res.append(vert_minus[i])
                else:
                    res.append([])
        return res, ind
    
    def solve_contact(self, points, i):
        for f in i:
            print(len(points[f]))
            m = sum(x.mass for x in points[f])
            point = sum((x.mass * x.pos) for x in points[f]) / m
            displacement = point - self.body.centre_of_mas_pos
            vel = self.body.velocity + np.cross(self.body.angular_velocity, displacement)
            relative_vel = np.dot(self.plane[f].unit_normal, vel) 
            denominator = np.cross(displacement, self.plane[f].unit_normal)
            denominator = np.matmul(denominator, np.linalg.inv(self.body.inertia_tensor).transpose())
            denominator = np.cross(denominator, displacement)
            denominator = np.dot(self.plane[f].unit_normal, denominator)
            denominator += 1 / self.body.mass
            j = -2 * relative_vel / denominator
            force = j * self.plane[f].unit_normal
            self.body.linear_momentum += force
            self.body.angular_momentum += np.cross(displacement, force)
            self.body.velocity = self.body.linear_momentum / self.body.mass
            self.body.angular_velocity = np.matmul(self.body.angular_momentum,
                                                   np.linalg.inv(self.body.inertia_tensor).transpose())

    def start_simulation(self, dt, fps):
        while True:
            time.sleep(1 / fps)
            #self.body.eilers_step(dt)
            self.body.midpoint_step(dt)
            c_p, number = self.collision_points()
            if number != []:
                self.solve_contact(c_p, number)
            self.body_view.update()


if __name__ == "__main__":
    cube = RigidBody(cube_generator(10, 0.2, [-10., 0., 30.]), [1., 50., 25., 100.], [-50., -50., 5000.], [0., 0., 25.])
    plane = []
    plane.append(Plane([[25, 25, -8], [25, -25, 1], [-25, -25, 10], [-25, 25, 1]]))
    plane.append(Plane([[60, 40, 60], [60, 40, -60], [-60, 40,-60], [-60, 40, 60]]))
    plane.append(Plane([[40, 60, 60], [40, 60, -60], [40, -60,-60], [40, -60, 60]]))
    plane.append(Plane([[60, -40, 60], [60, -40, -60], [-60, -40,-60], [-60, -40, 60]]))
    plane.append(Plane([[-40, 60, 60], [-40, 60, -60], [-40, -60,-60], [-40, -60, 60]]))
    s = Simulator(cube, plane)
