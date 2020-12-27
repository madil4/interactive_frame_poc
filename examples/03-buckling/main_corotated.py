import numpy as np
from readGTS import readGTS
from viewer import Viewer
from scipy.linalg import polar

# load mesh
x,t = readGTS("column2.gts")

# init
particles = np.shape(x)[0]
elements = np.shape(t)[0]
sc = np.zeros(elements) # for strain color
press = None;

# integration stuff
dt = 0.033
xi = x+0;
v = np.zeros((particles,2))
M = np.identity(particles * 2)
f_ex = np.zeros(particles*2)
for i in range(particles * 2):
    if i%2 == 0:
        f_ex[i] = 0
        f_ex[i+1] = -50


# # material 
u = 400

# boundary
cons = np.array([])
free = np.setdiff1d(np.arange(particles * 2), cons)


def per_draw(t):

	# time stepping
	xt = x.flatten()
	vt = v.flatten()

	xn = xt + 0
	y = xt + vt * dt

	# implicit
	# a = M[free].T[free].T - compute_k()[free].T[free].T * dt**2
	# b = M.dot(y)[free] + (compute_f(xt)[free] - compute_k()[free].T[free].T.dot(xt[free]) + f_ex[free]) * dt**2
	# xt[free] = np.linalg.inv(a).dot(b)

	# explicit
	xt[free] = y[free] + (compute_f(xt)[free]+f_ex[free]- 2.5*vt[free])*dt**2

	vt = (xt - xn) / dt

	xt = xt.reshape((particles,2))
	vt = vt.reshape((particles,2))

	for i in range(particles):

		if i == press:
			continue

		x[i] = xt[i]
		v[i] = vt[i]

		# ground collision
		if (x[i][1] < -10):
			x[i][1] = -10.0
			v[i][1] = -10.0


# ----- forces and stiffness
def compute_f(x):

	f = np.zeros(particles * 2)
	for i in range(elements):

		Dm = np.zeros((2, 2))
		Dm[:, 0] = xi[int(t[i][0])] - xi[int(t[i][2])]
		Dm[:, 1] = xi[int(t[i][1])] - xi[int(t[i][2])]

		s = np.zeros((2,particles))
		s[:, int(t[i][0])] = [1, 0]
		s[:, int(t[i][1])] = [0, 1]
		s[:, int(t[i][2])] = [-1, -1]

		I = np.identity(2)
		Bm = np.linalg.inv(Dm)
		g = np.kron(Bm.T.dot(s), I)

		w = abs(0.5 * np.linalg.det(Dm))

		# corotated
		F = g.dot(x).reshape(2,2)
		R,S = polar(F)
		P = np.zeros(4)
		F = F.reshape(4,)
		R= R.reshape(4,)
		P[0] = F[0] - R[0]
		P[1] = F[1] - R[1]
		P[2] = F[2] - R[2]
		P[3] = F[3]- R[3]
		s = abs(np.trace(S.reshape((2,2)) - I.reshape((2,2))))*10
		f_s = -w * g.T.dot(u*P)

		# # svtk
		# F = g.dot(x)
		# P = np.zeros(4)
		# P[0] = F[0]**3+F[0]*F[2]**2+F[0]*F[1]**2+F[2]*F[3]*F[1]-F[0]
		# P[1] = (F[0]**2)*F[1]+F[0]*F[2]*F[3]+F[1]**3+F[3]**2*F[1]-F[1]
		# P[2] = F[2]*F[0]**2+F[2]**3+F[3]*F[0]*F[1]+F[2]*F[3]**2-F[2]
		# P[3] = F[0]*F[1]*F[2]+(F[2]**2)*F[3]+(F[1]**2)*F[3]+F[3]**3-F[3]
		# s = abs(np.trace(F.reshape((2,2)) - I.reshape((2,2))))*10
		# f_s = -w * g.T.dot(u*P)

		f += f_s

		# for color mapping
		# s = abs(np.trace(F.reshape((2,2)) - I.reshape((2,2))))*10
		# s = abs(np.trace(S.reshape((2,2)) - I.reshape((2,2))))*10
		if s < 2:
			sc[i] = 1
		elif s > 2 and s < 10:
			sc[i] = 1-s/8
		else:
			sc[i] = 0

	return f

# ----- interactive

def on_press(event):
	global press
	
	l_min = 100;
	for i,x_p in enumerate(x):
		l = np.linalg.norm(x_p-[event.xdata,event.ydata])

		if l < 3:
			if l < l_min:
				l_min = l
				press = i

def on_motion(event):
	global press

	if press is None: return
	x[press] = [event.xdata,event.ydata]
	# x[press][1] = event.ydata

def on_release(event):
	global press

	press = None

Viewer.set_mesh(x,t,sc)
Viewer.set_per_draw(per_draw)
Viewer.set_events(on_press,on_release,on_motion)

Viewer.show()