

## Description


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for F+H<sub>2</sub>, and the 
corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For F+H<sub>2</sub>, we choose hyperspherical coordinate system, where &rho; is fixed at 10 Bohr. We employ other two 
coordinates, &theta; and &phi; to compute *ab initio* adiabatic PESs for three lowest electronic states and NACTs within them. The adiabatic 
potential energies are calculated in multireference configuration interaction (MRCI) level implementing Dunning’s correlation consistent 
aug-cc-pVTZ basis set with seven (7) electrons distributed over eight (8) active orbitals (7e, 8o) configuration active space (CAS). 
NACTs are computed using coupled-perturbed multi configuration space self-consistent field (CP-MCSCF) implemented in MOLPRO quantum chemistry
software. The first coordinate, &theta; ranges from 0 to &pi;/2 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number
of grid points will be 31 X 121. Integration is performed along 'Path 6' and plain text file format is used for both input and output.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_6'


>`adt num -mfile Tau_Theta.dat -nfile Tau_Phi.dat -efile Adiabatic_PES.dat -intpath 6 -txt`

---
