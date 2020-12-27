import numpy as np
from readGTS import readGTS
from viewer import Viewer

# load mesh
x,t = readGTS("beam.gts")

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
        f_ex[i+1] = 0

# material 
u = 527

# boundary
cons = np.array([0, 1, 2, 3])
free = np.setdiff1d(np.arange(particles * 2), cons)


def per_draw(t):

	# time stepping
	xt = x.flatten()
	vt = v.flatten()

	xn = xt + 0
	y = xt + vt * dt

	xt[free] = y[free] + (compute_f(xt)[free]+f_ex[free])*dt**2

	vt = (xt - xn) / dt

	xt = xt.reshape((particles,2))
	vt = vt.reshape((particles,2))

	for i in range(particles):

		if i == press:
			continue

		x[i] = xt[i]
		v[i] = vt[i]


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

		F = g.dot(x)

		e = np.array([F[0]-1,0.5*(F[1]+F[2]) ,0.5*(F[1]+F[2]) ,F[3]-1])		
		f_s = -w * g.T.dot(2*u * e)

		f += f_s

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
	# x[press] = [event.xdata,event.ydata]
	x[press][1] = event.ydata

def on_release(event):
	global press

	press = None

Viewer.set_mesh(x,t,sc)
Viewer.set_per_draw(per_draw)
Viewer.set_events(on_press,on_release,on_motion)

Viewer.show()