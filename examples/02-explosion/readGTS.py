import numpy as np

def readGTS(fname):

    fin = open(fname)

    # read the statistics
    linelist = fin.readline().split()

    V_size = int(linelist[0])
    F_size = int(linelist[1])

    # read the vertices
    V = np.zeros((V_size,2))
    for i in range(V_size):
        linelist = fin.readline().split()
        V[i] = np.array([float(linelist[0]), float(linelist[1])])

    # read the triangles
    F = np.zeros((F_size,3))
    for i in range(F_size):
        linelist = fin.readline().split()
        F[i] = np.array([int(linelist[0]),int(linelist[1]),int(linelist[2])])

    return V,F
