# Full Linear Elasticity

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

        F = g.dot(x).reshape(2, 2)
        I = np.identity(2)
        f_s = -w * g.T.dot((500 * (F + F.T - 2 * I) + 500 / 2 * np.trace(F - I) * I).reshape(4,))

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

        I_T = np.array([[2, 0., 0., 0],
                        [0, 1., 1., 0],
                        [0, 1., 1., 0],
                        [0, 0., 0., 2]])
        I_trans = np.array([[1, 0., 0., 1],
                            [0, 0., 0., 0],
                            [0, 0., 0., 0],
                            [1, 0., 0., 1]])
        k_s = -w * g.T.dot((500 * I_T + 500 / 2 * I_trans).dot(g))

        k += k_s

    return k


# Corotated Linear Elasticity

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

        k_s = -w * g.T.dot((2 * 500 * np.identity(4)).dot(g))

        k += k_s

    return k


# Simplest material model

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

        F = g.dot(x).reshape(2, 2)
        I = np.identity(2)
        f_s = -w * g.T.dot((500 * (F - I)).reshape(4,))

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

        k_s = -w * g.T.dot((500 * np.identity(4)).dot(g))

        k += k_s

    return k


# minimum surfaces
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

        F = g.dot(x).reshape(2, 2)
        f_s = -w * g.T.dot(g.dot(x))

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

        k_s = -w * g.T.dot((np.identity(4)).dot(g))

        k += k_s

    return k
