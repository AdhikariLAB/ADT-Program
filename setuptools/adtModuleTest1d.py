# An example to show use of ADT package through an python module

import numpy as np
from adt.numeric import adt_numeric


nact = np.loadtxt('../H3+Jacobi_test/tau_phi_mod.dat')
enr = np.loadtxt('../H3+Jacobi_test/energy_mod.dat')


grid  = nact[:,0]
ngrid = grid.shape
nstate = enr.shape[1]-1
ntau   = nact.shape[1]-1


angle, amat, db = adt_numeric.adt_quantities1d(grid, nact[:,1:], enr[:,1:])

angle = np.column_stack([grid, angle])
np.savetxt( 'angle.dat', angle ,delimiter="\t", fmt="%.8f")

for i in range(nstate):
    np.savetxt( 'diab_row_%s.dat'%(i+1), np.column_stack([grid, db[:,i,:]]) ,delimiter="\t", fmt="%.8f")

for i in range(nstate):
    np.savetxt( 'diab_row_%s.dat'%(i+1), np.column_stack([grid, amat[:,i,:]]) ,delimiter="\t", fmt="%.8f")