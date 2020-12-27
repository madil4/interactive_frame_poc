import numpy as np
from load_mesh import load_mesh
from view_mesh import view_mesh
from scipy.linalg import polar

# load mesh
mesh = load_mesh("beam.gts")


# needed function

def compute_f(x):

    f = np.zeros(mesh.ndof)
    for face in mesh.faces:

        Dm = np.zeros((2, 2))
        Dm[:, 0] = xi[face.vertex1.dof] - xi[face.vertex3.dof]
        Dm[:, 1] = xi[face.vertex2.dof] - xi[face.vertex3.dof]

        s = np.zeros((2, int(mesh.ndof / 2)))
        s[:, int(face.vertex1.dof[0] / 2)] = [1, 0]
        s[:, int(face.vertex2.dof[0] / 2)] = [0, 1]
        s[:, int(face.vertex3.dof[0] / 2)] = [-1, -1]

        w = abs(0.5 * np.linalg.det(Dm))
        Bm = np.linalg.inv(Dm)
        g = np.kron(Bm.T.dot(s), np.identity(2))

        # corotated - femdofo way
        F = g.dot(x).reshape(2, 2)
        R, S = polar(F)
        f_s = -w * g.T.dot(2 * 500 * (F - R).reshape(4,))

        f += f_s

    return f

def compute_k(x):

    k = np.zeros((mesh.ndof, mesh.ndof))
    for face in mesh.faces:

        Dm = np.zeros((2, 2))
        Dm[:, 0] = xi[face.vertex1.dof] - xi[face.vertex3.dof]
        Dm[:, 1] = xi[face.vertex2.dof] - xi[face.vertex3.dof]

        s = np.zeros((2, int(mesh.ndof / 2)))
        s[:, int(face.vertex1.dof[0] / 2)] = [1, 0]
        s[:, int(face.vertex2.dof[0] / 2)] = [0, 1]
        s[:, int(face.vertex3.dof[0] / 2)] = [-1, -1]

        w = abs(0.5 * np.linalg.det(Dm))
        Bm = np.linalg.inv(Dm)
        g = np.kron(Bm.T.dot(s), np.identity(2))

        # corotated - Ladislav way
        k_s = -w * g.T.dot((2 * 500 * np.identity(4)).dot(g))

        k += k_s

    return k


# init
dt = 0.033
cons = np.array([0, 1, 8, 9])
free = np.setdiff1d(np.arange(mesh.ndof), cons)
x = np.array([vertex.x for vertex in mesh.vertices]).flatten()
v = np.zeros(mesh.ndof)
xi = np.copy(x)
f_ex = np.zeros(mesh.ndof)

f_ex[[7, 15]] = -200

# compute Mass matrix
M = np.zeros((mesh.ndof, mesh.ndof))
for face in mesh.faces:

    Dm = np.zeros((2, 2))
    Dm[:, 0] = xi[face.vertex1.dof] - xi[face.vertex3.dof]
    Dm[:, 1] = xi[face.vertex2.dof] - xi[face.vertex3.dof]

    a = 0.01 * abs(0.5 * np.linalg.det(Dm)) / 3
    dofs = np.array([face.vertex1.dof, face.vertex2.dof, face.vertex3.dof]).flatten()
    M[np.ix_(dofs, dofs)] += a * np.identity(6)

M_inv = np.linalg.inv(M)


# solve
t = 10
xx = np.copy(x)
for i in range(int(t / dt)):

    # init
    xn = np.copy(x)

    # Backward Euler
    y = x + v * dt + M_inv.dot(f_ex) * dt**2
    a = M - compute_k(x) * dt**2
    b = M.dot(y - x) + compute_f(x) * dt**2
    dx = np.linalg.solve(a[free].T[free].T, b[free])
    x[free] += dx

    # velocity update
    v = (x - xn) / dt

    # store for animation
    xx = np.vstack((xx, x))


# view
view_mesh(xx, mesh.faces, dt)
