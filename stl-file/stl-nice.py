import numpy as np

# global variables
resfile = 'res.txt'
thetasteps = 50
hydrosteps = 200


# function to get the normal vector of a triangle, which is given by the position vector of the vertices
def trinormal(v1, v2, v3):
    normal = np.cross((v3 - v1), (v2 - v1))
    return normal


# function to calculate the surface area of a triangle based upon the normal vector
def areatri(nnn):
    area = np.linalg.norm(nnn) / 2
    return area


def writingsimplex(t1, t2, t3, nnn, fi):
    fi.write("facet normal " + str(nnn[0]) + " " + str(nnn[1]) + " "
             + str(nnn[2]) + "\n")
    fi.write("\touter loop \n")
    fi.write("\t\tvertex " + str(t1[0]) + " " + str(t1[1]) + " " + str(t1[2]) + " \n")
    fi.write("\t\tvertex " + str(t2[0]) + " " + str(t2[1]) + " " + str(t2[2]) + " \n")
    fi.write("\t\tvertex " + str(t3[0]) + " " + str(t3[1]) + " " + str(t3[2]) + " \n")
    fi.write("\tendloop \n")
    fi.write("endfacet\n")


a = np.loadtxt(resfile)
b = np.reshape(a, (hydrosteps, thetasteps , 3))

# print(b[0, 0, 0])
# print(b[0, 0, 1])
# print(b[0, 0, 2])

# print(b[0, 1, 0])
# print(b[0, 1, 1])
# print(b[0, 1, 2])

# print(b)
stlfile = resfile.split(".")[0] + ".stl"
f = open(stlfile, "w")
f.write("solid yield-surface \n")
currpos = np.zeros(3)
down1 = np.zeros(3)
down2 = np.zeros(3)
for i in range(1, hydrosteps):
    # posi = True
    for j in range(thetasteps):
        currpos[0:3] = b[i, j, 0:3]
        # if posi:
        # posi = False
        if j == thetasteps - 1:
            down1[0:3] = b[i - 1, 0, 0:3]
        else:
            down1[0:3] = b[i - 1, j + 1, 0:3]
        down2[0:3] = b[i - 1, j, 0:3]
        nn = trinormal(currpos, down1, down2)
        writingsimplex(t1=currpos, t2=down1, t3=down2, nnn=nn, fi=f)
        # else:
        # posi = True
    for j in range(thetasteps):
        currpos[0:3] = b[i, j, 0:3]
        down1[0:3] = b[i - 1, j, 0:3]
        if j == 0:
            down2[0:3] = b[i, thetasteps - 1, 0:3]
        else:
            down2[0:3] = b[i, j - 1, 0:3]
        nn = trinormal(currpos, down1, down2)
        writingsimplex(t1=currpos, t2=down1, t3=down2, nnn=nn, fi=f)

f.write("endsolid yield-surface \n")
f.close()
