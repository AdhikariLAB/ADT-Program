
## Description:


This folder contains five input files (molpro.config, atomfile.dat, geomfile.dat, frequency.dat and wilson.dat) of 
C<sub>6</sub>H<sub>3</sub>F<sub>3</sub>, and the output files corresponding to ADT quantities (ADT angles, residue of ADT angles, 
ADT matrices and diabatic potential energy matrices) calculated by this program package. For C<sub>6</sub>H<sub>3</sub>F<sub>3</sub> radical 
cation, we employ polar versions of two dimensionless normal mode coordinates, Q<sub>12x</sub> and Q<sub>12y</sub> (degenerate components of 
C-C-C in-plane scissoring vibration) to compute *ab initio* adiabatic potential energy surfaces (PESs) for six lowest electronic states and 
nonadiabatic coupling terms (NACTs) within them. The adiabatic potential energies are calculated in complete active space self-consistent 
field (CASSCF) level implementing Dunnings correlation consistent cc-pVDZ basis set and 11 electrons in 9 orbitals (11e, 9o) configuration 
active space (CAS). NACTs are computed using coupled-perturbed multi configuration space self-consistent field (CP-MCSCF) method implemented 
in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.02 to 0.06 and the second coordinate, &phi; spans from 0 to 
2&pi; such that total number of grid points will be 3 X 121. Integration is performed along 'Path 6' and HDF5 file format is used for output.

---
## Command:

The following command is used to generate the results, which are stored in the file 'ADT_numeric_irep_1_6.h5'


>`adt mol -config molpro.config -atomfile atomfile.dat -intpath 6 -h5`

---
