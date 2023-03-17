import vpython


def atov(arr):
    return vpython.vector(arr[0], arr[2], arr[1])


class RigidBodyView:
    def __init__(self, body):
        self.body = body
        self.mat_points = [vpython.sphere(pos=atov(particle.pos), radius=0.2)
                           for particle in body.particles]
        self.mat_points.append(
            vpython.sphere(pos=atov(self.body.centre_of_mas_pos), radius=0.2, color=vpython.color.red))

    def update(self):
        for point, particle in zip(self.mat_points, self.body.particles):
            point.pos = atov(particle.pos)
        self.mat_points[len(self.mat_points) - 1].pos = atov(self.body.centre_of_mas_pos)


class PlaneView:
    def __init__(self, plane):
        self.plane = []
        for i in plane[:2]:
            self.plane.append(
                vpython.quad(vs=[vpython.vertex(pos=atov(vert)) for vert in i.vertices], texture=vpython.textures.rug))

        for i in plane[2:4]:
            self.plane.append(vpython.quad(vs=[vpython.vertex(pos=atov(vert)) for vert in i.vertices],
                                           texture=vpython.textures.stucco))

        for i in plane[4:5]:
            self.plane.append(vpython.quad(vs=[vpython.vertex(pos=atov(vert)) for vert in i.vertices],
                                           texture=vpython.textures.metal))
