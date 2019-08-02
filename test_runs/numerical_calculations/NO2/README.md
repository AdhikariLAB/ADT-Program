

## Brief Description


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for NO<sub>2</sub>, and the 
corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For NO<sub>2</sub> radical, we employ polar versions of two dimensionless normal mode coordinates, Q<sub>1</sub> and 
Q<sub>3</sub> (bending and asymmetric stretching vibration) to compute *ab initio* adiabatic PESs for two low-lying electronic states (A<sup>'</sup>
states in C<sub>s</sub> geometry) and NACT within them. The adiabatic potential energies are calculated in complete active space self-consistent
field (CASSCF) level implementing Dunning’s correlation consistent polarized valence triple zeta (cc-pVTZ) basis set and 17 electrons in 11 
orbitals (17e, 11o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi configuration space self-consistent 
field (CP-MCSCF) method implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.2 to 14.0 and the second 
coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points will be 140 X 361. Integration is performed along 'Path 6' and for 
input and output, plain text file format is used.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_6'


>`adt num -mfile Tau_Rho.dat -nfile Tau_Phi.dat -efile Adiabatic_PES.dat -intpath 6`

---
