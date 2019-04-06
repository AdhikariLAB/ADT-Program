
## Description:


This folder contains the Adaiabtic PESs and NACTs for C<sub>6</sub>H<sub>3</sub>F<sub>3</sub> and the corresponding ADT Angles and Diabatic surfaces calculated by this program package.
For C<sub>6</sub>H<sub>3</sub>F<sub>3</sub> radical cation, we employ polar versions of two normal mode coordinates,
 Q<sub>10x</sub> and Q<sub>10y</sub> (degenerate components of 
C-C symmetric stretching vibration) to compute ab initio adiabatic potential energy surfaces (PESs) for six lowest electronic 
states and nonadiabatic coupling terms (NACTs) within them. The adiabatic potential energies are calculated in complete active 
space self-consistent field (CASSCF) level implementing Dunnings correlation consistent cc-pVDZ basis set and 11 electrons 
in 9 orbitals (11e, 9o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi configuration space 
self-consistent field (CP-MCSCF) method implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges 
from 0.0 to 1.0 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points wil be 48 X 180. 
Integration is performed along 'Path 6' and HDF5 file format is used for both input and output.

---
## Command:

The following command is used to generate the the results and the results are stored in the file 'ADT_numeric_6.h5'


>`adt num -nfile1 Tau_Rho.h5 -nfile2 Tau_Phi.h5 -efile Adiabatic_PES.h5 -intpath 6`

---