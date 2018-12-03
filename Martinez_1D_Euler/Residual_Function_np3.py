import numpy
import math

np = 3 # order of scheme

# Store the spectral coefficients for np = 3 and a fifth order polynomial:
volcoef = [0.0]*(np+1)

volcoef[0] = 1.0/12.0
volcoef[1] = 5.0/12.0
volcoef[2] = 5.0/12.0
volcoef[3] = 1.0/12.0

ncell = 10 #number of cells
nface = ncell+1 #cell boundaries
rk3 = 3 #This is because we are solving 3 equations(mass, momentum, energy)

Qf = numpy.zeros((rk3,np+1,ncell), dtype=float)

Qfl = numpy.empty((rk3),dtype=float)
Qfr = numpy.empty((rk3),dtype=float)

Fl = numpy.empty((rk3),dtype=float)
bflux = numpy.zeros((rk3, ncell, 2), dtype=float)

volflux = numpy.zeros((rk3,np+1,ncell), dtype=float)
h = numpy.zeros((rk3,np,ncell), dtype=float)

#f2c = numpy.empty((nface, 2), dtype=float)

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#---------------------------------------------------------------

def get_Qf(ncell, rk3, np, fle_gauss, Q):

	#Here we approximate U(x,t) = sum (basis_function(x)* U(t))
	
	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				Qf[iv,0,i] += fle_gauss[k,0]*Q[iv,k,i]
				Qf[iv,1,i] += fle_gauss[k,2]*Q[iv,k,i]
				Qf[iv,2,i] += fle_gauss[k,3]*Q[iv,k,i]
				Qf[iv,3,i] += fle_gauss[k,1]*Q[iv,k,i]
				
	return Qf
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#-----------BOUNDARY CONDITIONS---------------------------------

def get_Qbright_Qbleft(rk3,np, ncell, Qf):
	
	#Apply Boundary Conditions
	
	# Try to change the boundary conditions to see the impact to the solver
	#currently in comment block are the reflecting boundary conditions
	

	Qbleft = [0.0]*(rk3)
	Qbright = [0.0]*(rk3)
# 
# 	Qbleft[0] = Qf[0,0,0]
#  	Qbleft[1] = -Qf[1,0,0]
#  	Qbleft[2] = Qf[2,0,0]
#  
#  	Qbright[0] = Qf[0,np, ncell-1]
#  	Qbright[1] = -Qf[1, np, ncell-1]
#  	Qbright[2] = Qf[2, np, ncell-1]
	
	Qbleft[0] = 1.25
 	Qbleft[1] = 0.0
 	Qbleft[2] = 2.5
 
 	Qbright[0] = 0.125
 	Qbright[1] = 0.0
 	Qbright[2] = 0.25
	
	return Qbleft, Qbright
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#------------------BC FOR 1st CELL LEFT EDGE---------------------

def get_Qfl_Qfr_left(rk3,Qf,np, Qbleft):


	for iv in range(rk3):
		Qfl[iv] = Qbleft[iv]
		Qfr[iv] = Qf[iv,0,0]

	return Qfl, Qfr
#---------------------------------------------------------------

#Fl = get_rusanov_flux(Qfl, Qfr, normal, rk3)

def get_bflux_left(Fl):
	bflux[:,0,0] = Fl[:]
	return bflux
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#------------------BC FOR LAST CELL RIGHT EDGE-------------------

def get_Qfl_right(Qf,np, ncell):
	
	# Qfl = [1.0]*(rk3)
	
	Qfl[:] = Qf[:,np, ncell-1]
	return Qfl
	
#---------------------------------------------------------------
#Fl = get_rusanov_flux(Qfl, Qbright, normal, rk3)


def get_bflux_right(rk3, ncell, Fl):
	for i in range(rk3):
		bflux[i,ncell-1,1] = Fl[i]
	
	return bflux
	
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#------------------Interior Edge Flux----------------------------

def get_bflux_interior(nface, f2c, Qf, np, get_rusanov_flux):
	
	normal = 1.0
	
	ic_left = numpy.empty((nface-1), dtype=float)
	ic_right = numpy.empty((nface-1), dtype=float)
	Fl = numpy.empty((rk3), dtype = float)
	
	for i in range(1,nface-1):
		for j in range(rk3):
			
			ic_left[i] = f2c[i,0]
			ic_right[i] = f2c[i,1]
			
			ic_l = ic_left[i]
			ic_r = ic_right[i]
    
        	Qfl[j] = Qf[j, np, ic_l] 
        	Qfr[j] = Qf[j,0, ic_r]
		
		#----------------------------------#
		#Incorporate Rusanov flux function:
		
		"Returns the fluxes"
		gamma = 1.4
	
		rhol = numpy.abs(Qfl[0])
		ul   = Qfl[1]/rhol
		pl   = (gamma-1.0)*(Qfl[2]-(0.5*rhol*(ul**2)))
		
		rhor = numpy.abs(Qfr[0])
		ur   = Qfr[1]/rhor
		pr   = (gamma-1.0)*(Qfr[2]-(0.5*rhor*(ur**2)))

		
		c_a  = numpy.sqrt(((gamma*numpy.abs(pl+pr))/(rhol+rhor)))
		eigv = 0.5*(numpy.absolute(ul+ur)) + c_a 

		
		"Calculate the flux components"
	
		Fl[0] = Qfl[0]*ul + Qfr[0]*ur
		Fl[1] = Qfl[1]*ul + Qfr[1]*ur + (pl+pr)
		Fl[2] = (Qfl[2]+pl)*ul + (Qfr[2]+pr)*ur
		
		for iv in range(rk3):
			Fl[iv] = (Fl[iv]*normal - (Qfr[iv]-Qfl[iv]))*eigv*0.5
		
		#---------------------------------------#
		for j in range(rk3):
			bflux[j,ic_l,1]  = Fl[j]
			bflux[j,ic_r,0] = Fl[j]
    
	return bflux
	
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#----------------Compute All Volume Fluxes----------------------




def get_volflux(rk3, np, ncell, Qf):
	
	volflux = numpy.zeros((rk3,np+1,ncell), dtype=float)
	Fl = numpy.empty((rk3), dtype = float)
	
	for i in range(ncell):
		for j in range(np+1):
			for iv in range(rk3):
				Qfl[iv] = Qf[iv,j,i]
				
				gamma = 1.4
				rhol = numpy.abs(Qfl[0])
				ul   = Qfl[1]/rhol
				pl   = (gamma-1.0)*(Qfl[2]-(0.5*Qfl[0]*(ul**2)))
				
				Fl[0] = Qfl[0]*ul 
				Fl[1] = Qfl[1]*ul + pl
				Fl[2] = (Qfl[2]+pl)*ul
				
				
			for iv in range(rk3):
				volflux[iv,j,i] = Fl[iv]
	
	return volflux
        		

#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#-----------Get Initial h Matrix (RHS)--------------------------

def get_initial_hmatrix(rk3, np, ncell, bflux, fle_gauss):

	bl = [0.0]*(np)
	br = [0.0]*(np)

	for k in range(np):
		bl[k] = fle_gauss[k,0]
		br[k] = fle_gauss[k,1] 
		
	h = numpy.zeros((rk3,np,ncell), dtype=float)

	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				h[iv,k,i] = (bflux[iv,i,0]*bl[k]) - (bflux[iv,i,1]*br[k])
	return h
	
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------

#--------------Volume Integral by Gauss Lobatto----------------

def get_GL_hmatrix(rk3, np, ncell, h, volflux, volcoef, fle_grad_gauss):
	for i in range(ncell):
		for iv in range(rk3):
			for j in range(1, np):
				h[iv,j,i] = (h[iv,j,i] + 
                         	volflux[iv,0,i]*volcoef[0]*fle_grad_gauss[j,0] +
                         	volflux[iv,1,i]*volcoef[1]*fle_grad_gauss[j,2] +
                         	volflux[iv,2,i]*volcoef[2]*fle_grad_gauss[j,3] + 
                         	volflux[iv,3,i]*volcoef[3]*fle_grad_gauss[j,1] )
                         	
	return h 
		
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------	

#-------------Compute the Mass Matrix --------------------------

def get_residual_matrix(rk3, np, ncell, aic, h, dx):

	residual = numpy.zeros((rk3,np,ncell))
	
	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				residual[iv,k,i] = aic[k,k]*h[iv,k,i]/dx[i]
				
	return residual
	
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
		
def get_rusanov_flux(Ql, Qr, unitnorm, rk3):
	
	"Returns the fluxes"
	gamma = 1.4
	
	rhol = numpy.abs(Ql[0])
	ul   = Ql[1]/rhol
	pl   = (gamma-1.0)*(Ql[2]-(0.5*rhol*(ul**2)))
		
	rhor = numpy.abs(Qr[0])
	ur   = Qr[1]/rhor
	pr   = (gamma-1.0)*(Qr[2]-(0.5*rhor*(ur**2)))
		
	c_a  = numpy.sqrt(((gamma*numpy.abs(pl+pr))/(rhol+rhor)))
	eigv = 0.5*(numpy.absolute(ul+ur)) + c_a 

		
	"Calculate the flux components"
	
	Fl = numpy.empty((rk3), dtype = float)
	
	Fl[0] = Ql[0]*ul + Qr[0]*ur
	Fl[1] = Ql[1]*ul + Qr[1]*ur + (pl+pr)
	Fl[2] = (Ql[2]+pl)*ul + (Qr[2]+pr)*ur
		
	for iv in range(rk3):
		Fl[iv] = (Fl[iv]*unitnorm - (Qr[iv]-Ql[iv]))*eigv*0.5
	
	return Fl
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

#All volume flux:

def get_volume_flux(Ql):

	"Returns the volume fluxes"
	
	gamma = 1.4
	rhol = numpy.abs(Ql[0])
	ul   = Ql[1]/rhol
	pl   = (gamma-1.0)*(Ql[2]-(0.5*Ql[0]*(ul**2)))
	
	"Calculate the flux components"
	Fl[0] = numpy.abs(Ql[0]*ul) 
	Fl[1] = Ql[1]*ul + pl
	Fl[2] = (Ql[2]+pl)*ul 
	
	return Fl
#--------------------------------------------------------------------

