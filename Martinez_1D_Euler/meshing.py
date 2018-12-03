#Specify the number of faces for flux
import numpy
ncell = 10	
#------------------------------------
# set up the initial condition
xleft=  0.0
xlen= 1.0
xright= xlen + xleft
dxuni=xlen/ncell
dxmin=1.0
#------------------------------------
temp = numpy.empty((ncell+1), dtype=float)
for i in range(ncell+1):
    temp[i] = dxuni * (i)
       
#------------------------------------
def get_x_cell(ncell):
	x = numpy.empty((ncell), dtype=float)
	
	for i in range(ncell):
		x[i]=(temp[i+1]+(temp[i]))/2.0
	
	return x
	
#---------------------------------------

def get_dx_cell(ncell):
	dx = numpy.empty((ncell), dtype=float)
	
	for i in range(ncell):
		dx[i]=temp[i+1]-temp[i]
		
	return dx
	
#-----------------------------------------

def get_face_to_cell(nface):
	
	f2c = numpy.empty((nface, 2), dtype=float)
	for i in range(nface):
		f2c[i,0] = i-1
		f2c[i,1] = i
	f2c[nface-1,1] = 0.0
	f2c[0,0] = 0.0
	
	return f2c

