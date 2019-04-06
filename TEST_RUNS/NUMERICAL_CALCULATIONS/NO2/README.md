

## Brief Description


This folder contains the Adaiabtic PESs and NACTs for NO<sub>2</sub>, and the corresponding ADT Angles and Diabatic surfaces calculated by this program package. For NO<sub>2</sub> radical, we 
employ polar versions of two dimensionless normal mode coordinates, Q<sub>2</sub>  and Q<sub>3</sub>  (symmetric and asymmetric stretching 
vibrations) to compute ab initio adiabatic potential energy surfaces (PESs) for two lowest electronic states and nonadiabatic 
coupling terms (NACTs) within them. The adiabatic potential energies are calculated in multireference configuration 
interaction (MRCI) level implementing Dunning’s correlation consistent polarized valence triple zeta (cc-pVTZ) basis set and 
17 electrons in 10 orbitals (17e, 10o) configuration active space (CAS). NACTs are computed using numerical finite difference 
method (DDR) implemented in MOLPRO quantum chemistry software. The first coordinate, &rho; ranges from 0.013287 to 2.032944 
and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of grid points wil be 128 X 128. Integration
is performed along 'Path 1' and for input and output Plain Text file format is used.

---
## Command:

The following command is used to generate the the results and the results are stored in the folder 'ADT_numeric_1'


>`adt num -nfile1 Tau_Rho.dat -nfile2 Tau_Phi.dat -efile Adiabatic_PES.dat -intpath 1 -npy`

---