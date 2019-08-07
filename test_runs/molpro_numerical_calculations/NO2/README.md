

## Brief Description


This folder contains five input files (molpro.config, atomfile.dat, geomfile.dat, frequency.dat and wilson.dat) of NO<sub>2</sub>, and the
output files corresponding to ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) 
calculated by this program package. For NO<sub>2</sub> radical, we employ polar versions of two dimensionless normal mode coordinates, 
Q<sub>1</sub> and Q<sub>3</sub> (bending and asymmetric stretching vibrations) to compute *ab initio* adiabatic potential energy surfaces 
(PESs) for two lowest electronic states (A<sup>'</sup> states in C<sub>s</sub> geometry) and nonadiabatic coupling terms (NACTs) within them. 
The adiabatic potential energies are calculated in complete active space self-consistent field (CASSCF) level implementing Dunning’s 
correlation consistent polarized valence triple zeta (cc-pVTZ) basis set and 17 electrons in 11 orbitals (17e, 11o) configuration active 
space (CAS). NACTs are computed using coupled-perturbed multi-configuration space self-consistent field (CP-MCSCF) method implemented 
in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.2 to 12 and the second coordinate, &phi; spans from 0 to 
2&pi; such that total number of grid points will be 60 X 361. Integration is performed along 'Path 6' and plain text file format is used 
for output.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_irep_1_6'


>`adt mol -config molpro.config -atomfile atomfile.dat -intpath 6`

---
