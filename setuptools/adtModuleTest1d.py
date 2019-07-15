# An example to show use of ADT package as a python module

import numpy as np
from adt.numeric import adt_numeric as an
import os

nact = np.loadtxt('../H3+Jacobi_test/tau_phi_mod.dat')
enr = np.loadtxt('../H3+Jacobi_test/energy_mod.dat')


grid   = nact[:,0]
ngrid  = grid.shape
nstate = enr.shape[1]-1
ntau   = nact.shape[1]-1



fol = "ADT_numeric_script_1D"
os.makedirs(fol)

angle, res, amat, db = an.adt1d(grid, nact[:,1:], enr[:,1:])

angle = np.column_stack([grid, angle])
np.savetxt( '%s/Angle.dat'%fol, angle ,delimiter="\t", fmt="%.8f")


res = np.append(grid[-1], res)[None]  # <<<== to make it write as a row, as it just a single array
np.savetxt( '%s/Angle_residues.dat'%fol, res ,delimiter="\t", fmt="%.8f")

for i in range(nstate):
    np.savetxt( '%s/Diab_row_%s.dat'%(fol, i+1), np.column_stack([grid, db[:,i,:]]) ,delimiter="\t", fmt="%.8f")

for i in range(nstate):
    np.savetxt( '%s/Amat_row_%s.dat'%(fol, i+1), np.column_stack([grid, amat[:,i,:]]) ,delimiter="\t", fmt="%.8f")