import numpy as np


class Vertex():
    def __init__(self, x, y):
        self.x = np.array([x, y])
        self.dof = np.array([0, 0])


class Edge():
    def __init__(self, vertex1, vertex2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2


class Face():
    def __init__(self, edge1, edge2, edge3):
        self.edge1 = edge1
        self.edge2 = edge2
        self.edge3 = edge3

        self.vertex1 = edge1.vertex1
        self.vertex2 = edge1.vertex2

        # hot fix
        if edge2.vertex1.x[0] == edge1.vertex1.x[0] and edge2.vertex1.x[1] == edge1.vertex1.x[1]:
            self.vertex3 = edge2.vertex2
        elif edge2.vertex1.x[0] == edge1.vertex2.x[0] and edge2.vertex1.x[1] == edge1.vertex2.x[1]:
            self.vertex3 = edge2.vertex2
        else:
            self.vertex3 = edge2.vertex1


class Mesh():
    def __init__(self, vertices, edges, faces):
        self.vertices = vertices
        self.edges = edges
        self.faces = faces

        self.ndof = len(vertices) * 2


def load_mesh(fname):
    fin = open(fname)

    # read statistics
    linelist = fin.readline().split()

    nvertices = int(linelist[0])
    nedges = int(linelist[1])
    nfaces = int(linelist[2])

    # read vertices
    vertices = []
    j = 0
    for i in range(nvertices):
        linelist = fin.readline().split()

        node = Vertex(float(linelist[0]), float(linelist[1]))
        node.dof = [j, j + 1]
        j += 2

        vertices.append(node)

    # read edges
    edges = []
    for i in range(nedges):
        linelist = fin.readline().split()

        vertex1 = vertices[int(linelist[0]) - 1]
        vertex2 = vertices[int(linelist[1]) - 1]
        edge = Edge(vertex1, vertex2)

        edges.append(edge)

    # read faces
    faces = []
    for i in range(nfaces):
        linelist = fin.readline().split()

        edge1 = edges[int(linelist[0]) - 1]
        edge2 = edges[int(linelist[1]) - 1]
        edge3 = edges[int(linelist[2]) - 1]
        face = Face(edge1, edge2, edge3)

        faces.append(face)

    return Mesh(vertices, edges, faces)
