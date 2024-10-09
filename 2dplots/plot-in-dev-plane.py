import numpy as np
import matplotlib.pyplot as plt

pointonhydro = -100

filewithresults = '1res'

resopti = np.loadtxt(filewithresults)

saveresfile = 'results1.txt'

x = np.array([-1.0 / (2.0 ** 0.5), 1.0 / (2.0 ** 0.5), 0])

h = np.ones(3)

hn = h / np.linalg.norm(h)

xn = x / np.linalg.norm(x)

y = np.cross(a=hn, b=xn)

Psig = pointonhydro * hn

resres = np.zeros((len(resopti), 2))

count = 0

for i in resopti:

    r = i - Psig
    resres[count, 0] = np.dot(r, xn)
    resres[count, 1] = np.dot(r, y)

    count = count + 1

np.savetxt(fname=saveresfile, X=resres)
plt.plot(resres[:, 0], resres[:, 1], 'g-')
plt.savefig(filewithresults + '.png', dpi=600)
plt.show()
