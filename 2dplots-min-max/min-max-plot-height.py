import numpy as np
import matplotlib.pyplot as plt

filewithresults = 'res.txt'

saveresfile = 'data-height.txt'

thetasteps = 100

thetadir = 1

# min run
fff = filewithresults.split(".")[0] + "min." + filewithresults.split(".")[-1]

resopti = np.loadtxt(fff)

filelength = np.size(resopti[:, 0])

thetapoints = thetasteps 

hydropoints = int(filelength / thetapoints)

tdir1 = thetadir - 1

tdir2 = (tdir1 + int(thetapoints/2)) % thetapoints

if hydropoints != filelength / thetapoints:
    print("wrong steps in thetadir given")
    exit()

resres = np.zeros((hydropoints, 5))

sig1 = np.zeros(thetapoints)
sig2 = np.zeros(thetapoints)
sig3 = np.zeros(thetapoints)

h = np.ones(3)
hn = h / np.linalg.norm(h)

for i in range(hydropoints):
    lowerbound = i * thetapoints
    upperbound = (i + 1) * thetapoints
    
    # print(lowerbound)
    
    sig1 = resopti[lowerbound:upperbound, 0]
    sig2 = resopti[lowerbound:upperbound, 1]
    sig3 = resopti[lowerbound:upperbound, 2]
    
    x = np.sum(sig1) / np.size(sig1)
    y = np.sum(sig2) / np.size(sig2)
    z = np.sum(sig3) / np.size(sig3)
    
    hydpoint = (x + y + z) / 3 / hn[0]
    
    resres[i, 0] = hydpoint
    
    # print((i, np.size(sig1), hydpoint))
    
    Psig = hn * hydpoint
    p1 = resopti[lowerbound+tdir1,:]
    p2 = resopti[lowerbound+tdir2,:]
    
    r1 = p1 - Psig
    r2 = p2 - Psig
    
    resres[i, 1] = np.linalg.norm(r1)
    resres[i, 2] = -np.linalg.norm(r2)
    
maximilian =  np.max(np.absolute(resres[:, 1:3])) # max(np.max(resres[:, 1]), -np.min(resres[:, 2]))

resres[:, 3] = resres[:, 1] / maximilian
resres[:, 4] = resres[:, 2] / maximilian

sss = saveresfile.split(".")[0] + "min." + saveresfile.split(".")[-1]

np.savetxt(fname=sss, X=resres)

plt.plot(resres[:, 0], resres[:, 1], 'go-')
plt.plot(resres[:, 0], resres[:, 2], 'go-')
plt.savefig(filewithresults + 'min.png', dpi=600)
plt.show()


plt.plot(resres[:, 0], resres[:, 3], 'go-')
plt.plot(resres[:, 0], resres[:, 4], 'go-')
plt.savefig(filewithresults + '-relmin.png', dpi=600)
plt.show()

# max run
fff = filewithresults.split(".")[0] + "max." + filewithresults.split(".")[-1]

resopti = np.loadtxt(fff)

filelength = np.size(resopti[:, 0])

thetapoints = thetasteps 

hydropoints = int(filelength / thetapoints)

tdir1 = thetadir - 1

tdir2 = (tdir1 + int(thetapoints/2)) % thetapoints

if hydropoints != filelength / thetapoints:
    print("wrong steps in thetadir given")
    exit()

resres = np.zeros((hydropoints, 5))

sig1 = np.zeros(thetapoints)
sig2 = np.zeros(thetapoints)
sig3 = np.zeros(thetapoints)

h = np.ones(3)
hn = h / np.linalg.norm(h)

for i in range(hydropoints):
    lowerbound = i * thetapoints
    upperbound = (i + 1) * thetapoints
    
    # print(lowerbound)
    
    sig1 = resopti[lowerbound:upperbound, 0]
    sig2 = resopti[lowerbound:upperbound, 1]
    sig3 = resopti[lowerbound:upperbound, 2]
    
    x = np.sum(sig1) / np.size(sig1)
    y = np.sum(sig2) / np.size(sig2)
    z = np.sum(sig3) / np.size(sig3)
    
    hydpoint = (x + y + z) / 3 / hn[0]
    
    resres[i, 0] = hydpoint
    
    # print((i, np.size(sig1), hydpoint))
    
    Psig = hn * hydpoint
    p1 = resopti[lowerbound+tdir1,:]
    p2 = resopti[lowerbound+tdir2,:]
    
    r1 = p1 - Psig
    r2 = p2 - Psig
    
    resres[i, 1] = np.linalg.norm(r1)
    resres[i, 2] = -np.linalg.norm(r2)
    
maximilian =  np.max(np.absolute(resres[:, 1:3])) # max(np.max(resres[:, 1]), -np.min(resres[:, 2]))

resres[:, 3] = resres[:, 1] / maximilian
resres[:, 4] = resres[:, 2] / maximilian

sss = saveresfile.split(".")[0] + "max." + saveresfile.split(".")[-1]

np.savetxt(fname=sss, X=resres)

plt.plot(resres[:, 0], resres[:, 1], 'go-')
plt.plot(resres[:, 0], resres[:, 2], 'go-')
plt.savefig(filewithresults + 'max.png', dpi=600)
plt.show()


plt.plot(resres[:, 0], resres[:, 3], 'go-')
plt.plot(resres[:, 0], resres[:, 4], 'go-')
plt.savefig(filewithresults + '-relmax.png', dpi=600)
plt.show()

