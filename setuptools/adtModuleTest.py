# An example to show use of ADT package as a python module
import numpy as np
from adt.numeric import adt_numeric


def writeFile(data, file, xgrid):
    file = open(file,"wb")
    for tp in xgrid:
        np.savetxt( file, data[data[:,0]==tp] ,delimiter="\t", fmt="%.8f")
        file.write("\n")



nact1    = np.loadtxt('../TEST_RUNS/NUMERICAL_CALCULATIONS/H3+/Tau_Rho.dat')
nact2    = np.loadtxt('../TEST_RUNS/NUMERICAL_CALCULATIONS/H3+/Tau_Phi.dat')
enr      = np.loadtxt('../TEST_RUNS/NUMERICAL_CALCULATIONS/H3+/Adiabatic_PES.dat')
fullGrid = nact1[:,[0,1]]

grid1  = np.unique(nact1[:,0])
grid2  = np.unique(nact1[:,1])
ngrid1 = grid1.shape[0]
ngrid2 = grid2.shape[0]
nstate = enr.shape[1]-2
ntau   = nact1.shape[1]-2

nact1 = nact1[:,2:].reshape(ngrid1, ngrid2, ntau)
nact2 = nact2[:,2:].reshape(ngrid1, ngrid2, ntau)
enr   = enr[:,2:].reshape(ngrid1, ngrid2, nstate)


angle, res, amat,db = adt_numeric.adt_quantities(grid1, grid2, nact1, nact2, path = 1, energy = enr)


angle =  angle.reshape(ngrid1*ngrid2, ntau)
adtAngle = np.column_stack([fullGrid, angle])
writeFile(adtAngle, 'angle.dat', grid1)

res = np.column_stack([grid1, res])
np.savetxt( 'angle_residues.dat', res ,delimiter="\t", fmt="%.8f")

db = db.reshape(ngrid1*ngrid2, nstate, nstate)
for i in range(nstate):
    writeFile( np.column_stack([fullGrid, db[:,i,:]]), 'Diab_Row_%s.dat'%(i+1), grid1)

amat = amat.reshape(ngrid1*ngrid2, nstate, nstate)
for i in range(nstate):
    writeFile( np.column_stack([fullGrid, amat[:,i,:]]), 'Amat_Row_%s.dat'%(i+1), grid1)