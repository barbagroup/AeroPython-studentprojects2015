import numpy
from matplotlib import pyplot

fig, ax = pyplot.subplots(figsize=(40,2))

for t in range(1,200):
    if t == 1:
        fid = ''.join(['wake',str(t),'.dat'])
        with open(fid) as file_name:
            xw, yw= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
        points, = ax.plot(xw, yw, 'ro')
        fid = ''.join(['plate',str(t),'.dat'])
        with open(fid) as file_name:
            xp, yp= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
        coord, = ax.plot(xp, yp, 'b-', linewidth='2')

        #pyplot.grid(True)
        #pyplot.xlabel('x', fontsize=16)
        #pyplot.ylabel('y', fontsize=16)
        #pyplot.plot(xw, yw, 'ro')
        ax.set_xlim(-22,3)
        ax.set_ylim(-1,1)
    else:
        fid = ''.join(['wake',str(t),'.dat'])
        with open(fid) as file_name:
            new_x, new_y= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
        points.set_data(new_x, new_y)
        fid = ''.join(['plate',str(t),'.dat'])
        with open(fid) as file_name:
            new_x1, new_y1= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
        coord.set_data(new_x1, new_y1)
    pyplot.pause(0.1)

with open('wake200.dat') as file_name:
    xw, yw= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
points, = ax.plot(xw, yw, 'ro')
with open('plate200.dat') as file_name:
    xp, yp= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
coord, = pyplot.plot(xp, yp, 'b-', linewidth='2')
pyplot.show()
