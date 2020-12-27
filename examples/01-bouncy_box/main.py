import numpy as np
from readGTS import readGTS
from viewer import Viewer

# load mesh
x,t = readGTS("beam.gts")

# init
particles = np.shape(x)[0]
elements = np.shape(t)[0]
f_ex = np.zeros((particles,2))
v = np.zeros((particles,2))
xi = x+0;
sc = np.zeros(elements) # for strain color
press = None;


def per_draw(t):

	# compute gravity forces
	for i in range(particles):
		f_ex[i][0] = 0
		f_ex[i][1] = -50

	# compute elastic forces
	f = compute_f(x)

	for i in range(particles):

		if i == press:
			continue

		# time stepping
		v[i] += 0.033*(-f[i] + f_ex[i]) - 0.1*v[i]
		x[i] += 0.033*v[i]

		# ground collision
		if (x[i][1] < -10):
			x[i][1] = -10.0
			v[i][1] = -10.0


# ----- forces and stiffness
def compute_f(x):

	f = np.zeros(particles * 2)
	x = x.flatten()
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

		I = np.identity(2).flatten("F")
		w = abs(0.5 * np.linalg.det(Dm))

		F = g.dot(x)

		# for color mapping
		s = abs(np.trace(F.reshape((2,2)) - I.reshape((2,2))))*10
		if s < 2:
			sc[i] = 1
		elif s > 2 and s < 10:
			sc[i] = 1-s/8
		else:
			sc[i] = 0
		
		F_T = np.array([F[0], F[2], F[1], F[3]])

		F_I_trans = F[0] + F[3] - 2
		F_I_trans_I = np.array([F_I_trans, 0., 0., F_I_trans])

		f_s = w * g.T.dot(200 * (F + F_T - 2 * I) + 83.33 * F_I_trans_I)

		f += f_s

	return f.reshape(particles,2)


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

def on_release(event):
	global press

	press = None

Viewer.set_mesh(x,t,sc)
Viewer.set_per_draw(per_draw)
Viewer.set_events(on_press,on_release,on_motion)

Viewer.show()