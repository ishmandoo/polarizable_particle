import csv
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# imports data from Comsol
# the format is a whitespace separated table with x, y, z, E^2
# there is a data point every 0.05 um
# multiplying by 20 gives integers and allows me to store the field in a numpy array
lengthUnit = 0.05
lengthFactor = 1/lengthUnit

pointPos = np.array([5.,5.,2.])

eSquared = np.zeros((201,201,401))

def closestDistToPoint(point, line1, line2):
    return np.linalg.norm(np.cross((line2-line1),(line1-point)))/np.linalg.norm(line2-line1)

assert(closestDistToPoint(np.array([0,0,0]), np.array([1,0,1]), np.array([-1,0,1])) == 1)
assert(closestDistToPoint(np.array([0,0,0]), np.array([1,0,2]), np.array([-1,0,2])) == 2)
assert(closestDistToPoint(np.array([0,0,0]), np.array([1,0,0]), np.array([-1,0,0])) == 0)
assert(closestDistToPoint(np.array([0,0,0]), np.array([0,1,0]), np.array([0,1,1])) == 1)

class Particle:
    def __init__(self, x, y, z):
        self.polOverM = 6.11 * 10**-3
        self.p = np.array([x,y,z])
        self.v = np.array([0.,0.,-100000.])
        self.history = [self.p]
    def step(self, dt):
        # the acceleration grad(E^2) polarization / m should be in units of um/ms^2
        dE2 = getESquaredGrad(self.p)
        self.v += 0.5 * 0.5 * dE2 * self.polOverM * dt
        self.p += self.v * dt
        self.v += 0.5 * 0.5 * dE2 * self.polOverM * dt
        self.history.append(np.array(self.p))
    def steps(self, n, dt):
        for i in range(n):
            self.step(dt)


with open('esquared_fine_thin_pyramid.txt','rb') as tsvin:
    for row in tsvin:
        split = row.split()
        x = int(float(split[0])*lengthFactor)
        y = int(float(split[1])*lengthFactor)
        z = int(float(split[2])*lengthFactor)
        eSquared[x][y][z] = split[3]

def getESquared(pos):
    x = pos[0] * lengthFactor
    y = pos[1] * lengthFactor
    z = pos[2] * lengthFactor
    x1 = int(np.ceil(x))
    x0 = int(np.floor(y))
    y1 = int(np.ceil(y))
    y0 = int(np.floor(y))
    z1 = int(np.ceil(z))
    z0 = int(np.floor(z))

    print xp, x0, x1
    xd = (xp - x0)/(x1 - x0 + 0.001)
    yd = (yp - y0)/(y1 - y0 + 0.001)
    zd = (zp - z0)/(z1 - z0 + 0.001)

    c00 = (eSquared[x0][y0][z0] * (1 - xd)) + (eSquared[x1][y0][z0] * xd)
    c01 = (eSquared[x0][y0][z1] * (1 - xd)) + (eSquared[x1][y0][z1] * xd)
    c10 = (eSquared[x0][y1][z0] * (1 - xd)) + (eSquared[x1][y1][z0] * xd)
    c11 = (eSquared[x0][y1][z1] * (1 - xd)) + (eSquared[x1][y1][z1] * xd)

    c0 = (c00 * (1 - yd)) + (c10 * yd)
    c1 = (c01 * (1 - yd)) + (c11 * yd)

    c = (c0 * (1 - zd)) + (c1 * zd)
    return c

def getESquaredGrad(pos):
    # position units are in microns
    x = pos[0] * lengthFactor
    y = pos[1] * lengthFactor
    z = pos[2] * lengthFactor
    x0 = int(np.floor(x))
    y0 = int(np.floor(y))
    z0 = int(np.floor(z))

    #lengthUnitM = 10**6 * lengthUnit

    ddx = eSquared[x0 + 1][y0][z0] - eSquared[x0][y0][z0]
    ddy = eSquared[x0][y0 + 1][z0] - eSquared[x0][y0][z0]
    ddz = eSquared[x0][y0][z0 + 1] - eSquared[x0][y0][z0]
    return np.array([ddx/(lengthUnit * 1),ddy/(lengthUnit * 1),ddz/(lengthUnit * 1)])




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in np.linspace(4.9, 5.1, 5):
    #for j in np.linspace(4.9, 5.1, 5):
    part = Particle(i, 5, 2.3)
    print part.p
    while part.p[0] >= 0.0 and part.p[0] < 10.0 and part.p[1] >= 0.0 and part.p[1] < 10.0 and part.p[2] >= 2.0 and part.p[2] < 20.0:
        part.step(0.00000001)

    print closestDistToPoint(pointPos, part.history[-1], part.history[-2])
    trajectory = np.transpose(part.history)

    ax.scatter(trajectory[0],trajectory[1],trajectory[2])


'''
x = [4,6,5]
y = [5,5,5]
z = [1,1,2]
verts = [zip(x, y,z)]
ax.add_collection3d(Poly3DCollection(verts))

y = [4,6,5]
x = [5,5,5]
z = [1,1,2]
verts = [zip(x, y,z)]
ax.add_collection3d(Poly3DCollection(verts))
'''

ax.set_xlim([4.8,5.2])
ax.set_ylim([4.8,5.2])
ax.set_zlim([1,3])
plt.show()
