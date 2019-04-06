

## Description


This folder contains the Adaiabtic PESs and NACTs for H<sub>3</sub><sup>+</sup>, and the corresponding ADT Angles and Diabatic surfaces calculated by this program package. For H<sub>3</sub><sup>+</sup>, we 
choose hyperspherical coordinate system, where &rho; is fixed at 10 Bohr. We employ other two coordinates, &theta; and &phi; 
to compute ab initio adiabatic potential energy surfaces (PESs) for three lowest electronic states and nonadiabatic coupling 
terms (NACTs) within them. The adiabatic potential energies are calculated in multireference configuration interaction (MRCI) 
level implementing Dunning’s correlation consistent quintuple zeta (cc-pV5Z) basis set and 2 electrons in 10 orbitals (2e, 10o) 
configuration active space (CAS). NACTs are computed using numerical finite difference method (DDR) implemented in MOLPRO 
quantum chemistry software. The first coordinate, &theta; ranges from 0 to &pi;/2 and the second coordinate, &phi; spans 
from 0 to 2&pi; such that total number of grid points wil be 89 X 180. Integration is performed along 'Path 1' and Plain Text file format is used for both input and output

---
## Command:

The following command is used to generate the the results and the results are stored in the folder 'ADT_numeric_1'


>`adt num -nfile1 Tau_Rho.dat -nfile2 Tau_Phi.dat -efile Adiabatic_PES.dat -intpath 1 -txt`

---
