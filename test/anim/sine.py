import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation

x = np.linspace(-1,1,100)
y = np.sin(x)

fig = pl.figure()
ax = pl.axes(xlim=(-1,1),ylim=(-1,1))
line, = ax.plot([],[],lw=2)

def init():
    line.set_data(x,y)
    return line,

def animate(i):
    x = np.linspace(-1,1,1000)
    y = np.sin(2*np.pi*(x-0.01*i))
    line.set_data(x,y)
    return line,

anim = animation.FuncAnimation(fig,animate,init_func=init,frames=200,interval=20)

pl.show()

