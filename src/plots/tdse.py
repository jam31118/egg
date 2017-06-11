import matplotlib.pyplot as pl
from matplotlib import animation
import numpy as np
import parse

## Plotting
datafile = "./tdse.dat" #this should be modified to fetch this filename from CONFIG file
data = np.genfromtxt(datafile)

#print(type(len(data[0])))
#exit()

scalar = 1e2
for i in range(1,len(data)):
    data[i] = data[i]*scalar


print(len(data))
yNum = len(data) - 1
if (parse.tempmesh != yNum):
    print("yNum == "+str(yNum)+" tempmesh == "+str(parse.tempmesh))
    print("ERROR! not matching number of curves")
    exit(1)
dtau = parse.tmax / float(parse.tempmesh)
frameNum = parse.tempmesh 

fig = pl.figure()
ax = pl.axes(xlim=(parse.xmin,parse.xmax),
        ylim=(-0,1))

lines = []
line1, = ax.plot([],[],lw=1.5)
line2, = ax.plot([],[],':',lw=1.5)
lines.append(line1)
#lines.append(line2)

pl.xlabel(r"$\xi = \sqrt{m\omega/\hbar}x$",fontsize=14)
pl.ylabel(r"Amp.",fontsize=14)

time_text = ax.text(0.04,0.90,'',transform=ax.transAxes,fontsize=12)

def init():
    for line in lines:
        line.set_data([],[])
    time_text.set_text('')
    return tuple(lines) + (time_text,)

def animate(i):
    lines[0].set_data(data[0],data[i+1]**2)
    #y_sq = (data[i+1]) ** 2
    #lines[1].set_data(data[0],data[i+1]**2)
    tau = i*dtau
    time_text.set_text(r"$\tau$ =%5.2f" % tau) 
    return tuple(lines) + (time_text,)

anim = animation.FuncAnimation(fig,animate,init_func=init,
        frames=frameNum, interval=20, blit=False)

anim.save(parse.videoname + ".mp4",fps=30,extra_args=['-vcodec','libx264'])
#anim.save("doublewell_n=300.mp4",fps=30,extra_args=['-vcodec','libx264'])
#anim.save("doublewell_n=200.mp4",fps=30)

pl.show()


#for idx in range(yNum):
#    labelstr = "t = " + str(idx)
#    pl.plot(data[0],data[idx+1],label=labelstr)

#pl.legend()
#pl.show()

#pl.savefig("tdse.png") # should be modified to use parsed filename

