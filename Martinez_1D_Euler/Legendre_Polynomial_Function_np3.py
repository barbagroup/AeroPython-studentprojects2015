# Legendre Polynomial Function

import sympy
import numpy
from sympy.utilities.lambdify import lambdify

#-------------------------------------------------------
np = 3
xp = sympy.symbols('xp')
#-------------------------------------------------------

# Store the Quadrature Points and Weights for 4 point quadrature(np=3), fifth order polynomial:

#Quadrature Points

gauss = numpy.zeros((np+1, np+1), dtype=float)
gauss[0,0]=-0.5
gauss[1,0]=0.5
gauss[2,0]= -1 * numpy.sqrt(5)/10
gauss[3,0]= numpy.sqrt(5)/10.0

# Corresponding Quadrature weights

gauss[0,1]=1.0/12.0
gauss[1,1]=gauss[0,1]
gauss[2,1]=5.0/12.0
gauss[3,1]=gauss[2,1]

#----------------------------------------------------------
#Legendre Polynomial Terms
fle = [1.0]*(np+1)

fle[0] = 1.0
fle[1] = xp
fle[2] = (xp**2)-(1.0/12.0)
fle[3] = (xp**3)-(0.15*xp)
 
#---------------------------------------------------------
fle_eval = [1.0]*(np+1)

for i in range(np+1):
	fle_eval[i] = lambdify((xp),fle[i])
#----------------------------------------------------------
fle_0_eval = fle_eval[0]
fle_1_eval = fle_eval[1]
fle_2_eval = fle_eval[2]
fle_3_eval = fle_eval[3]

#----------------------------------------------------------
def fle0(xp):
	f0 = fle_0_eval(xp)
	return f0

def fle1(xp):
	f1 = fle_1_eval(xp)
	return f1

def fle2(xp):
	f2 = fle_2_eval(xp)
	return f2
    
def fle3(xp):
	f3 = fle_3_eval(xp)
	return f3

#----------------------------------------------------------------------	
#Legendre Polynomial Derivative
	
fle_grad = [0.0]*(np+1)
fle_grad[0] = 0.0
fle_grad[1] = 1
fle_grad[2] = 2*xp
fle_grad[3] = (3*(xp**2))-(0.15) 

#--------------------------------------------------------------------
fle_grad_eval = [1.0]*(np+1)

for i in range(np+1):
	fle_grad_eval[i] = lambdify((xp),fle_grad[i])
#--------------------------------------------------------------------
fle_grad_0_eval = fle_grad_eval[0]
fle_grad_1_eval = fle_grad_eval[1]
fle_grad_2_eval = fle_grad_eval[2]
fle_grad_3_eval = fle_grad_eval[3]


#--------------------------------------------------------------------
def fle0_grad(xp):
	fg0 = fle_grad_0_eval(xp)
	return fg0

def fle1_grad(xp):
	fg1 = fle_grad_1_eval(xp)
	return fg1

def fle2_grad(xp):
	fg2 = fle_grad_2_eval(xp)
	return fg2
    
def fle3_grad(xp):
	fg3 = fle_grad_3_eval(xp)
	return fg3

#----------------------------------------------------------------------
#-------------------------------------------------------------
def get_fle_gauss(np):
	"""
	fle_gauss => the Legendre Polynomial function evaluated at gauss points
			    	
	"""
	fle_ncb = [1.0]*(np+1)
	fle_ncb = [fle0(gauss[0,0]),fle1(gauss[0,0]),fle2(gauss[0,0]),fle3(gauss[0,0])]

	fle_pcb = [1.0]*(np+1)
	fle_pcb = [fle0(gauss[1,0]),fle1(gauss[1,0]),fle2(gauss[1,0]),fle3(gauss[1,0])]

	fle_nip = [1.0]*(np+1)
	fle_nip = [fle0(gauss[2,0]),fle1(gauss[2,0]),fle2(gauss[2,0]),fle3(gauss[2,0])]

	fle_pip = [1.0]*(np+1)
	fle_pip = [fle0(gauss[3,0]),fle1(gauss[3,0]),fle2(gauss[3,0]),fle3(gauss[3,0])]


#-------------------------------------------------------------------
	fle_gauss = numpy.zeros((np+1,np+1), dtype=float)

	for j in range(np+1):
		fle_gauss[j,0] = fle_ncb[j]
		fle_gauss[j,1] = fle_pcb[j]
		fle_gauss[j,2] = fle_nip[j]
		fle_gauss[j,3] = fle_pip[j]
		
	
	return fle_gauss
#-------------------------------------------------------------------
#-------------------------------------------------------------------
def get_fle_derivative(np):
	""" fle_grad_gauss => the Legendre Polynomial derivative of function evaluated
			   					 at gauss points 
	"""

	fle_grad_ncb = [1.0]*(np+1)
	fle_grad_ncb = [fle0_grad(gauss[0,0]),fle1_grad(gauss[0,0]),fle2_grad(gauss[0,0]),fle3_grad(gauss[0,0])]

	fle_grad_pcb = [1.0]*(np+1)
	fle_grad_pcb = [fle0_grad(gauss[1,0]),fle1_grad(gauss[1,0]),fle2_grad(gauss[1,0]),fle3_grad(gauss[1,0])]

	fle_grad_nip = [1.0]*(np+1)
	fle_grad_nip = [fle0_grad(gauss[2,0]),fle1_grad(gauss[2,0]),fle2_grad(gauss[2,0]),fle3_grad(gauss[2,0])]

	fle_grad_pip = [1.0]*(np+1)
	fle_grad_pip = [fle0_grad(gauss[3,0]),fle1_grad(gauss[3,0]),fle2_grad(gauss[3,0]),fle3_grad(gauss[3,0])]

#----------------------------------------------------------------------
	fle_grad_gauss = numpy.zeros((np+1,np+1), dtype=float)

	for j in range(np+1):
		fle_grad_gauss[j,0] = fle_grad_ncb[j]
		fle_grad_gauss[j,1] = fle_grad_pcb[j]
		fle_grad_gauss[j,2] = fle_grad_nip[j]
		fle_grad_gauss[j,3] = fle_grad_pip[j]
		
#-----------------------------------------------------------------------
	return fle_grad_gauss

  