

## Description


This folder contains two input files (molpro.config, atomfile.dat) of F+H<sub>2</sub>, and the output files corresponding to ADT quantities 
(ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this program package. 
For F+H<sub>2</sub>, we choose hyperspherical coordinate system, where &rho; is fixed at 5 Bohr. We employ other two 
coordinates, &theta; and &phi; to compute *ab initio* adiabatic potential energy surfaces (PESs) for three lowest electronic states and 
nonadiabatic coupling terms (NACTs) within them. The adiabatic potential energies are calculated in multireference configuration interaction 
(MRCI) level implementing Dunning’s correlation consistent aug-cc-pVTZ basis set with seven (7) electrons distributed over eight (8) active 
orbitals (7e, 8o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi configuration space self-consistent 
field (CP-MCSCF) implemented in MOLPRO quantum chemistry software. The first coordinate, &theta; ranges from 0 to &pi;/2 and the second coordinate, 
&phi; spans from 0 to 2&pi; such that total number of grid points will be 19 X 91. Integration is performed along 'Path 6' and HDF5 format is used 
for output.

---
## Command:

The following command is used to generate the results, which are stored in the file 'ADT_numeric_6.h5'


>`adt mol -config molpro.config -atomfile atomfile.dat -intpath 6 -h5`

---
