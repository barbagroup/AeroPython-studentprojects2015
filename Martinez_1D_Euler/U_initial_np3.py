import numpy

np = 3

aic = numpy.zeros((np, np), dtype=float)


# Store the spectral coefficients for np = 3 and a fifth order polynomial:
volcoef = [0.0]*(np+1)

volcoef[0] = 1.0/12.0
volcoef[1] = 5.0/12.0
volcoef[2] = 5.0/12.0
volcoef[3] = 1.0/12.0

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

#Inverse Mass Matrix
aic_3 = numpy.zeros((np,np), dtype=float)
aic_3[0,0]=1.0
aic_3[1,1]=12.0
aic_3[2,2]=180.0

#---------------------------------------
"""Defines inverse mass matrix and GL x points and weights
"""
def get_inverse_mass_matrix(np):
	
	#Third Order:Store the Inverse of the Mass Matrix
	
	for i in range(np):
		for j in range(np):
			aic[i,j] = aic_3[i,j]
	
	return aic
#------------------------------------------------------------------

def get_U_initial(np, rk3, ncell, fle_gauss, x , dx, aic):
	block = 0.50
	
	qinit = numpy.zeros((np+1,rk3,ncell), dtype=float)
	U_initial = numpy.empty((rk3,np,ncell), dtype=float)
	atemp = numpy.empty((np,rk3), dtype=float)
	
	for i in range(ncell):
		for j in range(np):
			for k in range(3):
			
				xeval = x[i] + gauss[j,0] * dx[i]
				if (xeval < block):
					qinit[k,0,i] = 1.25
					qinit[k,1,i] = 0.0
					qinit[k,2,i] = 2.5
				
				if (xeval >= block):
					qinit[k,0,i] = 0.125
					qinit[k,1,i] = 0.0
					qinit[k,2,i] = 0.25
					
				atemp[k,j] = ((qinit[0,k,i]*fle_gauss[j,0])*gauss[0,1] + (qinit[1,k,i]*fle_gauss[j,1])*gauss[1,1] +
                          	   (qinit[2,k,i]*fle_gauss[j,2])*gauss[2,1] + (qinit[3,k,i]*fle_gauss[j,3])*gauss[3,1])

			U_initial[0,j,i]=aic[j,j]*atemp[0,j]
			U_initial[1,j,i]=aic[j,j]*atemp[1,j]
			U_initial[2,j,i]=aic[j,j]*atemp[2,j]
				
	return U_initial
#-----------------------------------------------------------------------
