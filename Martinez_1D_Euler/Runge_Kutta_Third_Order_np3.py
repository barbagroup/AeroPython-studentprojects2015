import numpy          
#---------------------------------------------------
# ----------------Stage 1 Runge Kutta---------------

def get_Qrk_stage1(Q_rk, ncell, rk3, np, dt, residual):
	
	Qrk1 = numpy.zeros((rk3,np, ncell), dtype=float)
	
	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				Qrk1[iv,k,i] = Q_rk[iv,k,i]+(dt*residual[iv,k,i])
				
	return Qrk1
            
#------------------------------------------------------
#------------------------------------------------------
# ----------------Stage 1 Runge Kutta------------------

def get_Qrk_stage2(Q_rk, Qrk1, ncell, rk3, np, dt, residual):

	Qrk2 = numpy.zeros((rk3,np, ncell), dtype=float)
	
	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				Qrk2[iv,k,i] = (0.75*Q_rk[iv,k,i])+0.25*(Qrk1[iv,k,i]+(dt*residual[iv,k,i]))
				
	return Qrk2
            
#------------------------------------------------------
#------------------------------------------------------
# -----------------Stage 3 Runge Kutta-----------------

def get_Qrk_stage3(Q_rk, Qrk2, ncell, rk3, np, dt, residual):

	Qrk3 = numpy.zeros((rk3,np, ncell), dtype=float)
	
	for i in range(ncell):
		for iv in range(rk3):
			for k in range(np):
				Qrk3[iv,k,i] = ((1.0/3.0)*(Q_rk[iv,k,i])) + ((2.0/3.0)*(Qrk2[iv,k,i]+(dt*residual[iv,k,i])))
	
	return Qrk3
#------------------------------------------------------           
#ctime = ctime + 0.5*dt
#------------------------------------------------------
'''residnorm = 0.0

for i in range(ncell):
    for k in range(np):
        residnorm = residnorm + numpy.abs(residual[2,k,i])
        '''
