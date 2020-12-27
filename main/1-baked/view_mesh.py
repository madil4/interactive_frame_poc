from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.patches import Polygon


def view_mesh(x, faces, dt):

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 80), ylim=(0, 60))
    ax.grid()

    def loop(f):

        if f == 1:
            # init polygons
            for face in faces:
                polygon = Polygon([x[f][face.vertex1.dof], x[f][face.vertex2.dof], x[f][face.vertex3.dof]],
                                  zorder=2, edgecolor="C0", facecolor="white", linewidth=-1.5)
                ax.add_patch(polygon)

        elif f == len(x) - 1:
            # remove pologones
            [polygon.remove() for polygon in reversed(ax.patches)]

        elif f > 1:
            # update polygons
            for i, polygon in enumerate(ax.patches):
                polygon.set_xy([x[f][faces[i].vertex1.dof], x[f][faces[i].vertex2.dof], x[f][faces[i].vertex3.dof]])

    anim = animation.FuncAnimation(fig, loop, frames=len(x), interval=dt * 1000, repeat=True)
    plt.tight_layout(), plt.show()
