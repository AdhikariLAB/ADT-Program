

## Description:


This folder contains five input files (molpro.config, atomfile.dat, geomfile.dat, frequency.dat and wilson.dat) of NO<sub>3</sub>, and the 
output files corresponding to ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) 
calculated by this program package. For NO<sub>3</sub> radical, we employ polar versions of two dimensionless normal mode coordinates, 
Q<sub>4x</sub> and Q<sub>4y</sub> (degenerate components of in-plane asymmetric bending vibration) to compute *ab initio* adiabatic potential
energy surfaces (PESs) for five lowest electronic states and nonadiabatic coupling terms (NACTs) within them. The adiabatic potential 
energies are calculated in complete active space self-consistent field (CASSCF) level implementing 6-31G<sup>∗∗</sup> basis set and 9 
electrons in 8 orbitals (9e, 8o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi-configuration space 
self-consistent field (CP-MCSCF) method implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.1 to 
0.3 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points will be 3 X 121. Integration is performed 
along 'Path 6' and simple text format is used for output. 

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_irep_1_6'


>`adt mol -config molpro.config -atomfile atomfile.dat -intpath 6`

---
