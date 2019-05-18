# An example to show use of ADT package through an python module
import numpy as np
from adt.numeric import adt_numeric


nact1 = np.load('../TEST_RUNS/NUMERICAL_CALCULATIONS/NO3/Tau_Rho.npy')
nact2 = np.load('../TEST_RUNS/NUMERICAL_CALCULATIONS/NO3/Tau_Phi.npy')
enr = np.load('../TEST_RUNS/NUMERICAL_CALCULATIONS/NO3/Adiabatic_PES.npy')
fullGrid = nact1[:,[0,1]]

grid1  = np.unique(nact1[:,0])
grid2  = np.unique(nact1[:,1])
ngrid1 = 50
ngrid2 = 180
nstate = 5
ntau   = 10


nact1 = nact1[:,2:].reshape(ngrid1, ngrid2, ntau)
nact2 = nact2[:,2:].reshape(ngrid1, ngrid2, ntau)
enr = enr[:,2:].reshape(ngrid1, ngrid2, nstate)


angle, res, amat,db = adt_numeric.adt_quantities(grid1, grid2, nact1, nact2, path = 6, energy = enr)

print res.shape
print amat.shape
print db.shape
angle =  angle.reshape(ngrid1*ngrid2, ntau)

adtAngle = np.column_stack([fullGrid, angle])
# diab = diab.reshape(m*n,5,5)
# print diab[0]
np.savetxt( 'file', adtAngle ,delimiter="\t", fmt="%.8f")