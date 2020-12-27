import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.patches import Polygon


def view_mesh(x, faces, ex_loop, dt):

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 80), ylim=(0, 60))
    ax.grid()

    def loop(f):

        if f == 1:
            # init polygons
            for face in faces:
                polygon = Polygon([x[face.vertex1.dof], x[face.vertex2.dof], x[face.vertex3.dof]],
                                  zorder=2, edgecolor="C0", facecolor="white", linewidth=-1.5)
                ax.add_patch(polygon)

        elif f > 1:

            ex_loop(press, mouse)

            # update polygons
            for i, polygon in enumerate(ax.patches):
                polygon.set_xy([x[faces[i].vertex1.dof], x[faces[i].vertex2.dof], x[faces[i].vertex3.dof]])

    # for mouse evenet
    press = None
    mouse = None

    def on_press(event):
        nonlocal press, mouse

        l_min = 100
        for i, x_p in enumerate(x.reshape(-1, 2)):
            l = np.linalg.norm(x_p - [event.xdata, event.ydata])

            if l < 3:
                if l < l_min:
                    l_min = l
                    press = [2 * i, 2 * i + 1]

        mouse = [event.xdata, event.ydata]

    def on_motion(event):
        nonlocal press, x, mouse

        if press is None:
            return
        # x[press] = [event.xdata, event.ydata]
        mouse = [event.xdata, event.ydata]

    def on_release(event):
        nonlocal press

        press = None

    cid = fig.canvas.mpl_connect('button_press_event', on_press)
    cid = fig.canvas.mpl_connect('button_release_event', on_release)
    cid = fig.canvas.mpl_connect('motion_notify_event', on_motion)

    anim = animation.FuncAnimation(fig, loop, interval=dt * 1000)
    plt.tight_layout(), plt.show()
