import numpy as np

## Parsing CONFIG file
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
    elif (conf[0] == b"tempmesh"):
        tempmesh = int(conf[1].decode("utf-8")) + 1 # should be modified to easily understanable to human. Currently, it is not so clear when do I have to attach 1 or not.
    elif (conf[0] == b"xmax"):
        xmax = float(conf[1].decode("utf-8"))
        xmin = -xmax
    elif (conf[0] == b"tmax"):
        tmax = float(conf[1].decode("utf-8"))
    elif (conf[0] == b"videoname"):
        videoname = conf[1].decode("utf-8")

print("nodesStart:",nodesStart)
print("nodesEnd:",nodesEnd)


## Parsing energy.dat file
energyfile = "./energy.dat"
energydata = np.genfromtxt(energyfile,dtype="int,double",unpack=True)
#print(energydata[0][1])

