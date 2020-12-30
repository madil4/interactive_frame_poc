import numpy as np
from load_mesh import load_mesh
from view_mesh import view_mesh
from scipy.linalg import polar

# load mesh
mesh = load_mesh("frame.gts")


# needed function

def compute_r(x, g):

    r = np.zeros(4 * len(mesh.faces))
    for i, face in enumerate(mesh.faces):

        F = g[4 * i:4 * i + 4].dot(x).reshape(2, 2)
        R, S = polar(F)
        r_f = R.reshape(4,)

        r[4 * i:4 * i + 4] = r_f

    return r


# init
dt = 0.033
cons = np.array([0,1,2,3,46,47,48,49])
free = np.setdiff1d(np.arange(mesh.ndof), cons)
x = np.array([vertex.x for vertex in mesh.vertices]).flatten()
v = np.zeros(mesh.ndof)
xi = np.copy(x)
f_ex = np.zeros(mesh.ndof)

# compute g & k & M
g = np.zeros((4 * len(mesh.faces), mesh.ndof))
k = np.zeros((4 * len(mesh.faces), 4 * len(mesh.faces)))
M = np.zeros((mesh.ndof, mesh.ndof))
for i, face in enumerate(mesh.faces):

    Dm = np.zeros((2, 2))
    Dm[:, 0] = xi[face.vertex1.dof] - xi[face.vertex3.dof]
    Dm[:, 1] = xi[face.vertex2.dof] - xi[face.vertex3.dof]

    # compute g
    s = np.zeros((2, int(mesh.ndof / 2)))
    s[:, int(face.vertex1.dof[0] / 2)] = [1, 0]
    s[:, int(face.vertex2.dof[0] / 2)] = [0, 1]
    s[:, int(face.vertex3.dof[0] / 2)] = [-1, -1]

    Bm = np.linalg.inv(Dm)
    g_f = np.kron(Bm.T.dot(s), np.identity(2))

    g[4 * i:4 * i + 4] = g_f

    # compute k
    w = abs(0.5 * np.linalg.det(Dm))
    k_f = 2 * 1000 * dt**2 * w * np.identity(4)

    k[4 * i:4 * i + 4, 4 * i:4 * i + 4] = k_f

    # compute m
    # a = 0.0001* abs(0.5 * np.linalg.det(Dm)) / 3
    # dofs = np.array([face.vertex1.dof, face.vertex2.dof, face.vertex3.dof]).flatten()
    # M[np.ix_(dofs, dofs)] += a * np.identity(6)

M = np.identity(mesh.ndof)*1
M_inv = np.linalg.inv(M)


# solve
def ex_loop(press, mouse):
    global x, v, dt, free, cons,f_ex

    # freeze pressed
    if press is not None:
        # free = np.setdiff1d(free, press)
        f_ex[press] = 500 * (mouse - x[press])

    # init
    xn = np.copy(x)

    # Backward Euler
    a = M + g.T.dot(k.dot(g))
    y = x + v * dt + M_inv.dot(f_ex) * dt**2
    b = M.dot(y) + g.T.dot(k.dot(compute_r(x, g)))
    x[free] = np.linalg.inv(a).dot(b)[free]

    # velocity update
    v = (x - xn) / dt

    # free pressed
    if press is None:
        # free = np.setdiff1d(np.arange(len(x)), cons)
        f_ex = f_ex * 0 + 0


# view
view_mesh(x, mesh.faces, ex_loop, dt)
