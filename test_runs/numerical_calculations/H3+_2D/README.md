

## Description


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for H<sub>3</sub><sup>+</sup>, and 
the corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For H<sub>3</sub><sup>+</sup>, we choose hyperspherical coordinate system, where &rho; is fixed at 10 Bohr. We employ other 
two coordinates, &theta; and &phi; to compute *ab initio* adiabatic PESs for three lowest electronic states and NACTs within them. 
The adiabatic potential energies are calculated in multireference configuration interaction (MRCI) level implementing Dunning’s correlation 
consistent quintuple zeta (cc-pV5Z) basis set and 2 electrons in 10 orbitals (2e, 10o) configuration active space (CAS). NACTs are computed 
using numerical finite difference method (DDR) implemented in MOLPRO quantum chemistry software. The first coordinate, &theta; ranges from 
0 to &pi;/2 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points wil be 89 X 180. Integration is 
performed along 'Path 1' and Numpy binary format is used for both input and output.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_1'


>`adt num -mfile Tau_Theta.npy -nfile Tau_Phi.npy -efile Adiabatic_PES.npy -intpath 1 -nb`

---
