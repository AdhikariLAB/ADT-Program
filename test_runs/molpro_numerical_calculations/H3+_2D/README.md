

## Description


This folder contains two input files (molpro.config, atomfile.dat) of H<sub>3</sub><sup>+</sup>, and the output files corresponding to
ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For H<sub>3</sub><sup>+</sup>, we choose hyperspherical coordinate system, where &rho; is fixed at 4.5 Bohr. We employ other 
two coordinates, &theta; and &phi; to compute *ab initio* adiabatic potential energy surfaces (PESs) for three lowest electronic states and 
nonadiabatic coupling terms (NACTs) within them. The adiabatic potential energies are calculated in multireference configuration interaction 
(MRCI) level implementing Dunning’s correlation consistent quintuple zeta (cc-pV5Z) basis set and 2 electrons in 11 orbitals (2e, 11o) 
configuration active space (CAS). NACTs are computed using numerical finite difference method (DDR) implemented in MOLPRO quantum chemistry 
software. The first coordinate, &theta; ranges from &pi;/60 (3<sup>0</sup>) to &pi;/2 (90<sup>0</sup>) and the second coordinate, &phi; spans
from 0 to 2&pi; such that total number of grid points wil be 30 X 121. Integration is performed along 'Path 1' and Numpy binary format is used for 
both input and output.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_irep_1_1'


>`adt mol -config molpro.config -atomfile atomfile.dat -intpath 1 -nb`

---
