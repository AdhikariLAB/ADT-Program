import numpy as np
from adt.numeric import adt_numeric as an

nact = np.loadtxt('Tau_Phi.dat')
enr = np.loadtxt('Adiabatic_PES.dat')

grid = nact[:,0]
ngrid = grid.shape
nstate = enr.shape[1]-1
ntau = nact.shape[1]-1

angle, res, amat, db = an.adt1d(grid, nact[:,1:], enr[:,1:])
