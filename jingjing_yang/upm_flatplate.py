#import numpy
from matplotlib import pyplot
from functions_jy import *

rho = 1.225                                #kg/m^3
c   = 1.                                   #cord length
dx  = 0.01

np  = int(1/dx)                            #number of panels
x   = numpy.linspace(0.,1.,np+1)           #x-coord
y   = numpy.zeros(x.size,dtype=float)      #y-coord

T   = 20.                                  #simulation duration
nt  = 201                                  #number of time steps
dt  = T/(nt-1)                             #time step size
t   = numpy.linspace(0,T,nt)

#assume tailing edge is at origin when simulation begins
x0    = -1.*t                              #motion trajectory in x-dir
dx0dt = -1.                                #x velocity
y0    = numpy.zeros(x0.size,dtype=float)
dy0dt = 0.
theta = -5.*numpy.pi/180.
dthdt = 0.
gamma = numpy.zeros((np,nt-1),dtype=float) #vortex strength

#list of arrays for wake at different time step
W = [numpy.zeros((i,3),dtype=float) for i in range(1,nt)]

F = [numpy.zeros((np,3),dtype=float),W[0]]
#vortex locate at 1/4-length of each panel, V.M. Falker, 46'
for i in range(np):
    F[0][i,0] = 0.75*x[i] + 0.25*x[i+1]
    F[0][i,1] = 0.75*y[i] + 0.25*y[i+1]

#calculate the postion of first shedded vortex in the wake
[xw, yw] = first_vortex_pos(x,y,theta,dx0dt,dy0dt,dthdt,dt)

A    = coeff_matrix(x,y,xw,yw)
b    = rhs0(x,y,theta,dx0dt,dy0dt,dthdt)
solu = numpy.linalg.solve(A,b)

#transfer data for storage
F[0][0:np,2]  = solu[0:np]
F[1][0,0]       = xw
F[1][0,1]       = yw
F[1][0,2]       = solu[np]
gamma[0:np,0] = solu[0:np]

xylocal = numpy.empty((2,1),dtype=float)
xylocal[0:1,0] = F[1][0,0:1]
[xglob,yglob]= coord_trans(xylocal,x0[1],y0[1],theta,0)

W[0][0,0] = xglob
W[0][0,1] = yglob
W[0][0,2] = F[1][0,2]

#iteration
for i in range(2,nt):
    print('iter= ',i)
    j = i-1
    #xw, yw are constant if motion prescribed is
    #translational with constant speed
    [xw, yw] = first_vortex_pos(x,y,theta,dx0dt,dy0dt,dthdt,dt)

    #update F[1] to the current time step
    #basically store W[i-1] to F[1], then F[1]
    #F[1] will assume W[i-1]'s shape
    F[1] = numpy.empty((i,3),dtype=float)
    F[1][0,0] = xw; F[1][0,1] = yw;

    #xyglob assume W[i-2][:,0:2].T, #[:,inclusive:exclusive]
    xyglob = numpy.empty((2,i-1),dtype=float)
    xyglob[:,:] = W[i-2][:,0:2].T
    [xlocal,ylocal]= coord_trans(xyglob,x0[i],y0[i],theta,1)
    #if i==2: print(xlocal,ylocal)

    F[1][1:i,0] = xlocal; F[1][1:i,1] = ylocal; F[1][1:i,2] = W[i-2][:,2];
    #if i==2: print(F[1][1:i,:])

    A       = coeff_matrix(x,y,xw,yw)
    #total gamma of previous time step
    gamma_t = numpy.sum(gamma[:,i-2])
    #if i==2: print(gamma_t)

    b = rhs(x,y,theta,dx0dt,dy0dt,dthdt,F,gamma_t)
    #if i==2: print(b)

    solu = numpy.linalg.solve(A,b)
    F[0][0:np,2]  = solu[0:np]
    F[1][0,2]     = solu[np]
    #if i==2: print(F[0][0:np,2],F[1][0,2])

    F = iteration(x,y,theta,dx0dt,dy0dt,dthdt,F,dt,gamma_t)
    gamma[0:np,i-1] = F[0][0:np,2]

    #converting coordinates local -> global
    xylocal = F[1][:,0:2].T
    [xglob,yglob]= coord_trans(xylocal,x0[i],y0[i],theta,0)
    W[i-1][:,0] = xglob
    W[i-1][:,1] = yglob
    W[i-1][:,2] = F[1][:,2]

#end of iterations
#post-processing
animate_wake()

pyplot.figure(1,figsize=(10,0.5))
with open('./wake_profile/wake200.dat') as file_name:
    xw, yw= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
pyplot.plot(xw, yw, 'ro')
with open('./wake_profile/plate200.dat') as file_name:
    xp, yp= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
pyplot.plot(xp, yp, 'b-', linewidth='2')
pyplot.xlabel('x',fontsize=16)
pyplot.xlabel('y',fontsize=16)
pyplot.legend(['wake strucutre','flat plate'], loc='best')
pyplot.xlim(-23,3)
pyplot.ylim(-1,1)
pyplot.show()
