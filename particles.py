import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm as cm
import csv

class Field:
    def __init__(self, file, lengthUnit, origin, nx, ny, nz):
        self.lengthUnit = lengthUnit
        self.pointPos = np.array([0.5,0.5,2.])
        self.pointPosBot = np.array([0.5,0.5,1.95])
        self.origin = np.array(origin)
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.eSquared = self.readTsv(file, lengthUnit, nx, ny, nz)

    def indicesFromPosition(self, pos):
        return np.round((pos-self.origin)/self.lengthUnit).astype(int)

    def indicesAbovePosition(self, pos):
        return np.ceil((pos-self.origin)/self.lengthUnit).astype(int)

    def indicesBelowPosition(self, pos):
        return np.floor((pos-self.origin)/self.lengthUnit).astype(int)

    def positionFromIndex(self, index):
        return (index * self.lengthUnit) + self.origin

    def readCsv(self, file, lengthUnit, nx, ny, nz):
        eSquared = np.zeros((nx + 1,ny + 1,nz + 1))
        with open(file,'rb') as csvin:
            reader = csv.reader(csvin)
            for row in reader:
                x, y, z = self.indicesFromPosition(np.array(row[:3], dtype=float))
                eSquared[x][y][z] = row[3]
        return eSquared

    # imports data from Comsol
    # the format is a whitespace separated table with x, y, z, E^2
    # there is a data point every 0.05 um
    # multiplying by 20 gives integers and allows me to store the field in a numpy array
    def readTsv(self, file, lengthUnit, nx, ny, nz):
        eSquared = np.zeros((nx + 1,ny + 1,nz + 1))
        with open(file,'rb') as tsvin:
            for row in tsvin:
                split = row.split()
                x, y, z = self.indicesFromPosition(np.array(split[:3], dtype=float))
                assert(eSquared[x][y][z] == 0)
                eSquared[x][y][z] = split[3]
        return eSquared


    def plotField(self, y):
        fieldSlice = np.zeros((self.nx + 1,self.nz + 1))
        for i in range(self.nx):
            for j in range(self.nz):
                fieldSlice[i][j] = self.eSquared[i][self.indicesFromPosition(np.array([0,y,0]))[1]][j]

        plt.imshow(fieldSlice, cmap='hot', interpolation='nearest')
        plt.show()

    def nanSafe(self, num):
        if not np.isnan(num):
            return num
        else:
            return 0


    def getE(self, pos):
        x, y, z = pos / self.lengthUnit
        x1, y1, z1 = self.indicesAbovePosition(pos)
        x0, y0, z0 = self.indicesBelowPosition(pos)

        xd = (x - x0)/(x1 - x0 + 0.001)
        yd = (y - y0)/(y1 - y0 + 0.001)
        zd = (z - z0)/(z1 - z0 + 0.001)

        if not all([not np.isnan(e2) for e2 in [self.eSquared[x0][y0][z0],self.eSquared[x1][y0][z0],self.eSquared[x0][y0][z1],self.eSquared[x1][y0][z1],self.eSquared[x0][y1][z0],self.eSquared[x1][y1][z0],self.eSquared[x0][y1][z1],self.eSquared[x1][y1][z1]]]):
            return 0

        c00 = (self.eSquared[x0][y0][z0] * (1 - xd)) + (self.eSquared[x1][y0][z0] * xd)
        c01 = (self.eSquared[x0][y0][z1] * (1 - xd)) + (self.eSquared[x1][y0][z1] * xd)
        c10 = (self.eSquared[x0][y1][z0] * (1 - xd)) + (self.eSquared[x1][y1][z0] * xd)
        c11 = (self.eSquared[x0][y1][z1] * (1 - xd)) + (self.eSquared[x1][y1][z1] * xd)

        c0 = (c00 * (1 - yd)) + (c10 * yd)
        c1 = (c01 * (1 - yd)) + (c11 * yd)

        c = (c0 * (1 - zd)) + (c1 * zd)
        return np.sqrt(c)

    # returns the gradient in the E^2 field at position pos
    def getESquaredGrad(self, pos):
        x0, y0, z0 = self.indicesBelowPosition(pos)

        # caluclates differences between pos and the next field value over in x, y, and z
        # returns vector of differences over step length
        # units of esquared are um^2*kg^2/(ms^6*A^2) so this function returns units of um*kg^2/(ms^6*A^2)
        ddx = self.eSquared[x0 + 1][y0][z0] - self.eSquared[x0][y0][z0]
        ddy = self.eSquared[x0][y0 + 1][z0] - self.eSquared[x0][y0][z0]
        ddz = self.eSquared[x0][y0][z0 + 1] - self.eSquared[x0][y0][z0]
        return np.array([ddx/self.lengthUnit, ddy/self.lengthUnit, ddz/self.lengthUnit])


class Particle:
    def __init__(self, field, p, v):
        self.field = field
        # polarization / m in units of ms^4 A^2 kg^-2
        self.polOverM = 6.11 * 10**-3
        self.p = p
        self.v = v
        self.history = [self.p]
        self.peakField = 0
    # advances the particle by a time dt
    def step(self, dt):

        # check if position is valid
        if all(np.logical_not(np.isnan(self.p))):
            self.peakField = max(self.peakField, self.field.getE(self.p))

        # the acceleration grad(E^2) polarization / m should be in units of um/ms^2
        dE2 = self.field.getESquaredGrad(self.p)

        print 0.5 * dE2 * self.polOverM * dt
        self.v += 0.5 * 0.5 * dE2 * self.polOverM * dt
        self.p += self.v * dt
        self.v += 0.5 * 0.5 * dE2 * self.polOverM * dt
        self.history.append(np.array(self.p))
    def steps(self, n, dt):
        for i in range(n):
            self.step(dt)
    def closestApproachSegment(self, line1, line2):
        point = self.field.pointPos
        return np.linalg.norm(np.cross((line2-line1),(line1-point)))/np.linalg.norm(line2-line1)
    def closestApproachLineLine(self, line11, line12, line21, line22):
        a = line12 - line11
        b = line22 - line21
        c = line21 - line11
        if (np.isnan(a).any()) or (np.isnan(b).any()) or (np.isnan(c).any()):
            return 1000
        #print a, b, c
        return np.linalg.norm(np.dot(c, np.cross(a, b)))/np.linalg.norm(np.cross(a, b))
    def closestApproach(self):
        #return np.nanmin(np.array([self.closestApproachLineLine(self.history[i], self.history[i - 1], self.field.pointPos, self.field.pointPosBot) for i in range(1,len(self.history))]))
        closest_index = np.nanargmin([np.linalg.norm(self.field.pointPos-pos) for pos in self.history])
        if closest_index >= len(self.history)-1:
            return self.closestApproachLineLine(self.history[closest_index - 1], self.history[closest_index], self.field.pointPos, self.field.pointPosBot)
        return np.fmin(self.closestApproachLineLine(self.history[closest_index], self.history[closest_index + 1], self.field.pointPos, self.field.pointPosBot), self.closestApproachLineLine(self.history[closest_index - 1], self.history[closest_index], self.field.pointPos, self.field.pointPosBot))


class ParticleSet:
    def __init__(self, field, ps, v, bounds):
        self.particles = [Particle(field, np.array(p), np.array(v)) for p in ps]
        self.xMin = bounds[0]
        self.xMax = bounds[1]
        self.yMin = bounds[2]
        self.yMax = bounds[3]
        self.zMin = bounds[4]
        self.zMax = bounds[5]
        for part in self.particles:
            print part.p
            print part.v
    def run(self, stepSize):
        for part in self.particles:
            while self.isInBounds(part):
                part.step(stepSize)
    def isInBounds(self, part):
        return part.p[0] >= self.xMin and part.p[0] < self.xMax and part.p[1] >= self.yMin and part.p[1] < self.yMax and part.p[2] >= self.zMin and part.p[2] < self.zMax
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')


        x = [.49,.51,.5]
        y = [.5,.5,.5]
        z = [1,1,2]
        verts = [zip(x, y,z)]
        ax.add_collection3d(Poly3DCollection(verts))


        colors = cm.rainbow(np.linspace(0, 1, len(self.particles)))
        for part, color in zip(self.particles, colors):
            trajectory = np.transpose(part.history)

            ax.scatter(trajectory[0],trajectory[1],trajectory[2], color=color)


        ax.set_xlim([self.xMin,self.xMax])
        ax.set_ylim([self.yMin,self.yMax])
        ax.set_zlim([self.zMin,self.zMax])
        plt.show()

    def printFinalPositions(self):
        for part in self.particles:
            print part.p

    def printPeakField(self):
        for part in self.particles:
            print part.peakField

    def printClosestApproaches(self):
        for part in self.particles:
            print part.closestApproach()



#field = Field("C:\\Users\\ishma\\Dropbox (SteinLab)\\simulation\\small_spacing_pyramid\\field.txt", 0.005, 200, 200, 600)

#parts = ParticleSet(field, [[x,x,2.5] for x in np.linspace(.49, .51, 5)], [0.,0.,-100000], [0., 1., 0., 1., 1., 3.])
#parts.run(0.000000001)
#parts.printFinalPositions()
#parts.printPeakField()
#parts.plot()
