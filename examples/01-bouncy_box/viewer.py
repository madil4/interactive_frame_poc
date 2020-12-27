import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib import cm
import time


class Viewer:

    cm = cm.get_cmap('hot', 50)

    @classmethod
    def set_mesh(cls,x,t,sc):
        cls.x = x
        cls.t = t
        cls.sc = sc

    @classmethod
    def set_per_draw(cls,per_draw):
        cls.per_draw = per_draw

    @classmethod
    def set_events(cls,on_press,on_release,on_motion):
        cls.on_press = on_press
        cls.on_release = on_release
        cls.on_motion = on_motion

    @classmethod
    def show(cls):

        # the loop function
        def framing(f):
            
            t1 = time.time()
            
            if f==1:

                # render polygon
                cls.polygons = []
                for f in cls.t:
                    polygon = Polygon([cls.x[int(f[0])],cls.x[int(f[1])],cls.x[int(f[2])]],zorder=2,edgecolor="C0",facecolor=cls.cm(0.0),linewidth=-1.5)
                    ax.add_patch(polygon)
                    
                    # store for updates
                    cls.polygons.append(polygon)
                    
            elif f>1:
                
                # set the loop
                cls.per_draw(f)

                # update polygon
                for i,polygon in enumerate(cls.polygons):
                    polygon.set_xy([ cls.x[int(cls.t[i][0])] ,
                                     cls.x[int(cls.t[i][1])] ,
                                     cls.x[int(cls.t[i][2])] ])

                    polygon.set_facecolor(cls.cm(cls.sc[i]))

            t2 = time.time()
            # print(1000*(t2 - t1))


        # view the mesh
        fig = plt.figure()
        ax = plt.axes(xlim=(-10, 40), ylim=(-25, 25), aspect="equal")

        # for mouse evenet
        cid = fig.canvas.mpl_connect('button_press_event', cls.on_press )
        cid = fig.canvas.mpl_connect('button_release_event', cls.on_release )
        cid = fig.canvas.mpl_connect('motion_notify_event', cls.on_motion )

        # the loop for debug use ,frames=3,repeat=False
        anim = animation.FuncAnimation(fig, framing,interval=33)


        line = Line2D((-10,40),(-10,-10))
        ax.add_line(line)

        plt.grid()
        plt.show()


