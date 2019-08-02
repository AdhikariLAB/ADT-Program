

## Description:


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for 
C<sub>6</sub>H<sub>6</sub><sup>+</sup>, and the corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic 
potential energy matrices) calculated by this program package. For C<sub>6</sub>H<sub>6</sub><sup>+</sup> radical cation, we 
employ polar versions of two dimensionless normal mode coordinates, Q<sub>16x</sub> and Q<sub>16y</sub> (degenerate components of in-plane 
asymmetric stretching vibration) to compute *ab initio* adiabatic PESs for five lowest electronic states and NACTs within them. The adiabatic
potential energies are calculated in complete active space self-consistent field (CASSCF) level implementing 6-31G<sup>∗∗</sup> basis set 
and 29 electrons in 15 orbitals (29e, 15o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi configuration 
space self-consistent field (CP-MCSCF) method implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.0 
to 2.0 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points wil be 20 X 180. Integration
is performed along 'Path 6' and Numpy Binary format is used for both input and output. 

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_6'


>`adt num -mfile Tau_Rho.npy -nfile Tau_Phi.npy -efile Adiabatic_PES.npy -intpath 6 -nb`

---
