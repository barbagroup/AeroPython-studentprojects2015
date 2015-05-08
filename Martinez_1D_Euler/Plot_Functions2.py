
# Legendre Polynomial Function for Plots

import sympy
import numpy
from sympy.utilities.lambdify import lambdify

#-------------------------------------------------------
np = 3
xp = sympy.symbols('xp')
#-------------------------------------------------------


#Quadrature Points

gauss_plot = numpy.zeros((np), dtype=float)
gauss_plot[0]=-1./3.
gauss_plot[1]=0.
gauss_plot[2]=1./3.

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


fle_np = [1.0]*(np+1)
fle_np = [fle0(gauss_plot[0]),fle1(gauss_plot[0]),fle2(gauss_plot[0]),fle3(gauss_plot[0])]

fle_zero = [1.0]*(np+1)
fle_zero = [fle0(gauss_plot[1]),fle1(gauss_plot[1]),fle2(gauss_plot[1]),fle3(gauss_plot[1])]

fle_pp = [1.0]*(np+1)
fle_pp = [fle0(gauss_plot[2]),fle1(gauss_plot[2]),fle2(gauss_plot[2]),fle3(gauss_plot[2])]

#-------------------------------------------------------------------
fle_gauss_plot = numpy.zeros((np+1,np), dtype=float)

for j in range(np+1):
	fle_gauss_plot[j,0] = fle_np[j]
	fle_gauss_plot[j,1] = fle_zero[j]
	fle_gauss_plot[j,2] = fle_pp[j]

#-----------------------------------------------------------------------------		
def plot_density2(ncell,np, Q, x, dx):
	
	
	rho1 = numpy.zeros((ncell), dtype=float)
	rho2 = numpy.zeros((ncell), dtype=float)
	rho3 = numpy.zeros((ncell), dtype=float)
	
	xx1 = numpy.zeros((ncell), dtype=float)
	xx2 = numpy.zeros((ncell), dtype=float)
	xx3 = numpy.zeros((ncell), dtype=float)
	
	for i in range(ncell):
		for j in range(np):
			
			rho1[i] += fle_gauss_plot[j,0]*Q[0,j,i]
			rho2[i] += fle_gauss_plot[j,1]*Q[0,j,i]
			rho3[i] += fle_gauss_plot[j,2]*Q[0,j,i]
			
		xx1[i] = x[i] + (-1./3.)*dx[i]
		xx2[i] = x[i] + (0.)*dx[i]
		xx3[i] = x[i] + (1./3.)*dx[i]

		
	rho_t2 = numpy.concatenate((rho1,rho2, rho3))
	xx_t2 = numpy.concatenate((xx1, xx2, xx3))
		
	return xx_t2, rho_t2
	
#-----------------------------------------------------------------------------		
def plot_velocity2(ncell,np, Q, x, dx):

	rho1 = numpy.zeros((ncell), dtype=float)
	rho2 = numpy.zeros((ncell), dtype=float)
	rho3 = numpy.zeros((ncell), dtype=float)
	
	mom1 = numpy.zeros((ncell), dtype=float)
	mom2 = numpy.zeros((ncell), dtype=float)
	mom3 = numpy.zeros((ncell), dtype=float)
	
	vel1 = numpy.zeros((ncell), dtype=float)
	vel2 = numpy.zeros((ncell), dtype=float)
	vel3 = numpy.zeros((ncell), dtype=float)
	
	xx1 = numpy.zeros((ncell), dtype=float)
	xx2 = numpy.zeros((ncell), dtype=float)
	xx3 = numpy.zeros((ncell), dtype=float)
	
	for i in range(ncell):
		for j in range(np):
			
			rho1[i] += fle_gauss_plot[j,0]*Q[0,j,i]
			rho2[i] += fle_gauss_plot[j,1]*Q[0,j,i]
			rho3[i] += fle_gauss_plot[j,2]*Q[0,j,i]
			
			mom1[i] += fle_gauss_plot[j,0]*Q[1,j,i]
			mom2[i] += fle_gauss_plot[j,1]*Q[1,j,i]
			mom3[i] += fle_gauss_plot[j,2]*Q[1,j,i]
		
		vel1[i] = mom1[i]/rho1[i]
		vel2[i] = mom2[i]/rho2[i]
		vel3[i] = mom3[i]/rho3[i]
			
		xx1[i] = x[i] + (-1./3.)*dx[i]
		xx2[i] = x[i] + (0.)*dx[i]
		xx3[i] = x[i] + (1./3.)*dx[i]
		
	vel_t2 = numpy.concatenate((vel1,vel2, vel3))
	xx_t2 = numpy.concatenate((xx1, xx2, xx3))
		
	return xx_t2, vel_t2
	
#------------------PLOT Pressure!-----------------------------------------------	

def plot_pressure2(ncell,np, Q, x, dx):
	
	gamma = 1.4
	
	rho1 = numpy.zeros((ncell), dtype=float)
	rho2 = numpy.zeros((ncell), dtype=float)
	rho3 = numpy.zeros((ncell), dtype=float)
	
	mom1 = numpy.zeros((ncell), dtype=float)
	mom2 = numpy.zeros((ncell), dtype=float)
	mom3 = numpy.zeros((ncell), dtype=float)
	
	en1 = numpy.zeros((ncell), dtype=float)
	en2 = numpy.zeros((ncell), dtype=float)
	en3 = numpy.zeros((ncell), dtype=float)
	
	press1 = numpy.zeros((ncell), dtype=float)
	press2 = numpy.zeros((ncell), dtype=float)
	press3 = numpy.zeros((ncell), dtype=float)
		
	xx1 = numpy.zeros((ncell), dtype=float)
	xx2 = numpy.zeros((ncell), dtype=float)
	xx3 = numpy.zeros((ncell), dtype=float)
	
	for i in range(ncell):
		for j in range(np):
			
			rho1[i] += fle_gauss_plot[j,0]*Q[0,j,i]
			rho2[i] += fle_gauss_plot[j,1]*Q[0,j,i]
			rho3[i] += fle_gauss_plot[j,2]*Q[0,j,i]
			
			mom1[i] += fle_gauss_plot[j,0]*Q[1,j,i]
			mom2[i] += fle_gauss_plot[j,1]*Q[1,j,i]
			mom3[i] += fle_gauss_plot[j,2]*Q[1,j,i]
			
			en1[i] += fle_gauss_plot[j,0]*Q[2,j,i]
			en2[i] += fle_gauss_plot[j,1]*Q[2,j,i]
			en3[i] += fle_gauss_plot[j,2]*Q[2,j,i]
			
		press1[i] = (gamma-1.)*(en1[i]- 0.5*(mom1[i]**2)/rho1[i]) 
		press2[i] = (gamma-1.)*(en2[i]- 0.5*(mom2[i]**2)/rho2[i]) 
		press3[i] = (gamma-1.)*(en3[i]- 0.5*(mom3[i]**2)/rho3[i]) 
			
		xx1[i] = x[i] + (-1./3.)*dx[i]
		xx2[i] = x[i] + (0.)*dx[i]
		xx3[i] = x[i] + (1./3.)*dx[i]
		
	press_t2 = numpy.concatenate((press1,press2, press3))
	xx_t2 = numpy.concatenate((xx1, xx2, xx3))
		
	return xx_t2, press_t2
	
#------------------PLOT Energy!-----------------------------------------------	

def plot_energy2(ncell,np, Q, x, dx):	

	en1 = numpy.zeros((ncell), dtype=float)
	en2 = numpy.zeros((ncell), dtype=float)
	en3 = numpy.zeros((ncell), dtype=float)
	
	xx1 = numpy.zeros((ncell), dtype=float)
	xx2 = numpy.zeros((ncell), dtype=float)
	xx3 = numpy.zeros((ncell), dtype=float)
	
	for i in range(ncell):
		for j in range(np):
			
			en1[i] += fle_gauss_plot[j,0]*Q[2,j,i]
			en2[i] += fle_gauss_plot[j,1]*Q[2,j,i]
			en3[i] += fle_gauss_plot[j,2]*Q[2,j,i]
	
		xx1[i] = x[i] + (-1./3.)*dx[i]
		xx2[i] = x[i] + (0.)*dx[i]
		xx3[i] = x[i] + (1./3.)*dx[i]
	
	en_t2 = numpy.concatenate((en1,en2, en3))
	xx_t2 = numpy.concatenate((xx1, xx2, xx3))
		
	return xx_t2, en_t2

	