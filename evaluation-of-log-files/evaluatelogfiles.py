import numpy as np

numberprocs = 5

logfilename = ""
logfileextension = ""

solver = "nr"  # es 

timeincs = []
totaltime = []
tries = []

for i in range(1, numberprocs + 1):
    a = logfilename + str(i) + logfileextension
    f = open(a, "r")
    lines = f.readlines()
    for j in lines:
        timeincs.append(float((j.split("time increment:")[-1]).split("total time:")[0]))
        totaltime.append(float(j.split("total time:")[-1]))
        if solver == "nr":
            tries.append(int((j.split("with:")[-1]).split("tries")[0]))

if tries:
    a = np.stack((tries, timeincs, totaltime), 1)
else:
    a = np.stack((timeincs, totaltime), 1)

np.savetxt("infologfile.dat", a)
