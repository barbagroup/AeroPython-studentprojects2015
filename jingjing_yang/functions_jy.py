import numpy
from matplotlib import pyplot

def first_vortex_pos(x,y,theta,u0,v0,dthdt,dt):
    """k, fraction for the distance of the 1st vortex
          to the trailing edge
       M, transformation matrix
    """
    k = 0.2

    r  = numpy.sqrt((x[-1]-x[-2])**2+(y[-1]-y[-2])**2)
    tx = (x[-1]-x[-2])/r
    ty = (y[-1]-y[-2])/r

    #get vel on the last panel at t=0
    M = numpy.empty((2,2))
    M = numpy.array( [ [ numpy.cos(theta), numpy.sin(theta)],
                       [-numpy.sin(theta), numpy.cos(theta)] ])

    r1 =-M.dot( numpy.array([u0,v0]) )
    r2 =-dthdt*numpy.array([-(y[-1]-0.),x[-1]-0.])
    #r1 = [0.9961947, 0.08715574]
    #r2 = [0., 0.]

    r  = r1+r2
    u1 = r[0]
    v1 = r[1]

    #Distance of the first vortex from the trailing edge
    dh=k*dt*(u1*tx+v1*ty)

    xw=x[-1]+dh*tx
    yw=y[-1]+dh*ty

    return [xw, yw]

def coeff_matrix(x,y,xw,yw):

    n = x.size
    np= n-1
    A = numpy.zeros((n,n),dtype=float)

    A[-1,:] = numpy.ones(n)

    for i in range(np):

        #control point at 3-quarter of each panel, V.M. Falker, 46'
        xn = 0.25*x[i]+0.75*x[i+1]
        yn = 0.25*y[i]+0.75*y[i+1]

        #normal of i-th panel
        r = numpy.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)

        nx= -(y[i+1]-y[i])/r
        ny=  (x[i+1]-x[i])/r

        for j in range(np):

            #vortex locate at 1/4-length of each panel, V.M. Falker, 46'
            xc = 0.75*x[j]+0.25*x[j+1]
            yc = 0.75*y[j]+0.25*y[j+1]

            #velocity induced on the control point of i-th panel by the j-th vortex
            [u,v] = vortex2d(xn,yn,xc,yc,1.)

            A[i,j] = u*nx + v*ny


        #influence in the wake
        [u,v] = vortex2d(xn,yn,xw,yw,1.)

        A[i,-1] = u*nx+v*ny

    return A

def vortex2d(x0,y0,x1,y1,gamma):

    r = (x0-x1)**2+(y0-y1)**2
    u=0.
    v=0.

    #use linear formula for samll distance
    delta = 1e-8
    if r > delta:
        u = 0.5*gamma*(y0-y1)/(numpy.pi*r);
        v =-0.5*gamma*(x0-x1)/(numpy.pi*r);
    else:
        u = 0.5*gamma*(y0-y1)/(numpy.pi*delta);
        v =-0.5*gamma*(x0-x1)/(numpy.pi*delta);

    return [u,v]

def rhs0(x,y,theta,u0,v0,dthdt):
    """
    create the rhs for the first time step
    """

    n = x.size
    np= n-1
    b=numpy.zeros(n);

    for i in range(np):
        #control point at 3-quarter of each panel, V.M. Falker, 46'
        xn=0.25*x[i]+0.75*x[i+1]
        yn=0.25*y[i]+0.75*y[i+1]

        r = numpy.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)
        nx =-(y[i+1]-y[i])/r
        ny = (x[i+1]-x[i])/r

        #get vel on the last panel at t=0
        M = numpy.empty((2,2))
        M = numpy.array( [ [ numpy.cos(theta), numpy.sin(theta)],
                           [-numpy.sin(theta), numpy.cos(theta)] ])

        r1 =-M.dot( numpy.array([u0,v0]) )
        r2 =-dthdt*numpy.array([-(yn-0.),xn-0.])
        r  = r1+r2

        u = r[0]
        v = r[1]

        b[i] =-u*nx-v*ny

    return b

def coord_trans(xy1,x0,y0,theta,m):

    M = numpy.empty((2,2))
    if m==1 :

        M = numpy.array( [ [ numpy.cos(theta), numpy.sin(theta)],
                           [-numpy.sin(theta), numpy.cos(theta)] ])
        xy1[0,:] -= x0
        xy1[1,:] -= y0
        r = M.dot(xy1)

    elif m==0 :

        M = numpy.array( [ [ numpy.cos(theta), -numpy.sin(theta)],
                           [ numpy.sin(theta),  numpy.cos(theta)] ])
        r1=M.dot(xy1)
        i, j=r1.shape
        r2 = numpy.empty((2,j),dtype=float)
        r2[0,:] = x0
        r2[1,:] = y0
        r=r1+r2

    x=r[0,:]; y=r[1,:]
    return [x,y]


def rhs(x,y,theta,u0,v0,dthdt,F,gamma_t):

    n      = x.size
    np     = n-1
    #[m,cc]=size(F{1,2});
    [m,cc] = F[1].shape
    b      = numpy.zeros(n)
    b[n-1] = gamma_t

    for i in range(np):

        xn=0.25*x[i]+0.75*x[i+1]
        yn=0.25*y[i]+0.75*y[i+1]

        r = numpy.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)

        nx=-(y[i+1]-y[i])/r; ny=(x[i+1]-x[i])/r;

        #get vel on the last panel at t=0
        M = numpy.empty((2,2))
        M = numpy.array( [ [ numpy.cos(theta), numpy.sin(theta)],
                           [-numpy.sin(theta), numpy.cos(theta)] ])

        r1 =-M.dot( numpy.array([u0,v0]) )
        r2 =-dthdt*numpy.array([-(yn-0.),xn-0.])
        r  = r1+r2

        u = r[0]
        v = r[1]


        u1=0; v1=0;
        for j in range(1,m):
            [u2,v2] = vortex2d(xn,yn,F[1][j,0],F[1][j,1],F[1][j,2])
            u1=u1+u2; v1=v1+v2;

        u = u+u1; v = v+v1;

        b[i]=-u*nx-v*ny;

    return b

def iteration(x,y,theta,u0,v0,dthdt,F,dt,gamma_t):

    k  = 0.33
    it = 3 #number of iteration

    [r,c]=F[1].shape

    #updating first vortex (1=yes, 0=no)
    pv=0;

    for m in range(it):
        #updating position of vortex in the wake
        for n in range(1-pv,r):

            #position of n-th vortex in the wake
            xc=F[1][n,0]; yc=F[1][n,1];

            #velocity induced by whirling system
            [u,v]=vel1(xc,yc,F);

            #updating
            F[1][n,0]=xc+k*dt*u; F[1][n,1]=yc+k*dt*v;

        #new solution
        #calculating position of the first vortex in the wake
        xw=F[1][0,0]; yw=F[1][0,1];

        A    = coeff_matrix(x,y,xw,yw)
        b    = rhs(x,y,theta,u0,v0,dthdt,F,gamma_t)
        solu = numpy.linalg.solve(A,b)

        n1 =solu.size
        np=n1-1

        F[0][0:np,2]  = solu[0:np]
        F[1][0,2]     = solu[np]

    return F

def vel1(x,y,F):

    u=0; v=0;
    for i in range(2):
        [m,n]=F[i].shape

        for j in range(m):
            xc   =F[i][j,0]
            yc   =F[i][j,1]
            gamma=F[i][j,2]

            [u1,v1] = vortex2d(x,y,xc,yc,gamma)
            u=u+u1; v=v+v1;

    return [u,v]

def animate_wake():
    fig, ax = pyplot.subplots(figsize=(40,2))

    for t in range(1,200):
        if t == 1:
            fid = ''.join(['./wake_profile/wake',str(t),'.dat'])
            with open(fid) as file_name:
                xw, yw= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
            points, = ax.plot(xw, yw, 'ro')
            fid = ''.join(['./wake_profile/plate',str(t),'.dat'])
            with open(fid) as file_name:
                xp, yp= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
            coord, = ax.plot(xp, yp, 'b-', linewidth='2')

            ax.set_xlim(-22,3)
            ax.set_ylim(-1,1)
        else:
            fid = ''.join(['./wake_profile/wake',str(t),'.dat'])
            with open(fid) as file_name:
                new_x, new_y= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
            points.set_data(new_x, new_y)
            fid = ''.join(['./wake_profile/plate',str(t),'.dat'])
            with open(fid) as file_name:
                new_x1, new_y1= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
            coord.set_data(new_x1, new_y1)
        pyplot.pause(0.1)

    with open('./wake_profile/wake200.dat') as file_name:
        xw, yw= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
    points, = ax.plot(xw, yw, 'ro')
    with open('./wake_profile/plate200.dat') as file_name:
        xp, yp= numpy.loadtxt(file_name, dtype=float, delimiter=',', unpack=True)
    coord, = pyplot.plot(xp, yp, 'b-', linewidth='2')
    pyplot.show()
