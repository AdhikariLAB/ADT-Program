import numpy as np
from adt.numeric import adt_numeric as an

nact1 = np.load('Tau_Rho.npy')
nact2 = np.load('Tau_Phi.npy')
enr = np.load('Adiabatic_PES.npy')

fullGrid = nact1[:,[0,1]]
grid1 = np.unique(nact1[:,0])
grid2 = np.unique(nact1[:,1])
ngrid1 = grid1.shape[0]
ngrid2 = grid2.shape[0]
nstate = enr.shape[1]-2
ntau = nact1.shape[1]-2

nact1 = nact1[:,2:].reshape(ngrid1, ngrid2, ntau)
nact2 = nact2[:,2:].reshape(ngrid1, ngrid2, ntau)
enr = enr[:,2:].reshape(ngrid1, ngrid2, nstate)

angle, res, amat, db = an.adt2d(grid1, grid2, nact1, nact2, path = 6, enr)
