import matplotlib.pyplot as pl
import numpy as np
import parse

"""
## Parse CONFIG
configfile = "./CONFIG"
config = np.genfromtxt(configfile,dtype="S20,S20",delimiter="=")

nodes = ""
nodesStart = -1
nodesEnd = -1
for conf in config:
    if (conf[0] == b"nodes"):
        nodes = conf[1].decode("utf-8")

nodesStart = int(nodes.split(":")[0])
nodesEnd = int(nodes.split(":")[1])
nodesNum = nodesEnd - nodesStart + 1

print("nodesStart:",nodesStart)
print("nodesEnd:",nodesEnd)
"""

## Plotting
datafile = "./y.dat"
#datafile = "./tdse.dat"
data = np.genfromtxt(datafile)

print(len(data))
yNum = len(data) - 1
if (parse.nodesNum is not yNum):
    print("ERROR! not matching nodesNum")
    exit()

for idx in range(yNum):
    print(parse.energydata[idx])
    energy = float(parse.energydata[idx][1])
    energyStr = "%.3e" % energy
    labelstr = "e = " + energyStr
    nodesNo = parse.nodesStart + idx
    pl.plot(data[0],data[nodesNo],label=labelstr)
pl.xlabel(r"$\xi = (mK/\hbar^2)^{(1/4)} x$",fontsize=14)
pl.ylabel(r"Amp.",fontsize=14)
pl.legend()
pl.show()

#pl.savefig("doublewell.png") # should be modified to use parsed filename
