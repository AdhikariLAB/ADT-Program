# An example to show use of ADT package through an python module

import numpy as np
from adt.numeric import adt_numeric


nact = np.loadtxt('/home/koushik/CODES/adt-program/H3+Jacobi_test/tau_phi_mod.dat')
enr = np.loadtxt('/home/koushik/CODES/adt-program/H3+Jacobi_test/energy_mod.dat')


grid  = nact[:,0]
ngrid = grid.shape
nstate = enr.shape[1]-1
ntau   = nact.shape[1]-1


angle, _, db = adt_numeric.adt_quantities1d(grid, nact[:,1:], enr[:,1:])

angle = np.column_stack([grid, angle])
np.savetxt( 'file', angle ,delimiter="\t", fmt="%.8f")