

## Description:


This folder contains the Adaiabtic PESs and NACTs for NO<sub>3</sub> and the corresponding ADT Angles and Diabatic surfaces calculated by this program package. For NO<sub>3</sub> radical, we 
employ polar versions of two normal mode coordinates, Q<sub>3x</sub> and Q<sub>3y</sub> (degenerate components of in-plane asymmetric 
stretching vibration) to compute ab initio adiabatic potential energy surfaces (PESs) for five lowest electronic states and 
nonadiabatic coupling terms (NACTs) within them. The adiabatic potential energies are calculated in complete active space 
self-consistent field (CASSCF) level implementing 6-31G<sup>∗∗</sup> basis set and 9 electrons in 8 orbitals (9e, 8o) configuration 
active space (CAS). NACTs are computed using coupled-perturbed multi configuration space self-consistent field (CP-MCSCF) 
method implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.0 to 5.0 and 
the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points wil be 50 X 180. Integration
is performed along 'Path 6' and Numpy Binary format is used for both input and ouput. 

---
## Command:

The following command is used to generate the the results and the results are stored in the folder 'ADT_numeric_6'


>`adt num -nfile1 Tau_Rho.npy -nfile2 Tau_Phi.npy -efile Adiabatic_PES.npy -intpath 6 -nb`

---