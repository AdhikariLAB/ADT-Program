

## Description:


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for NO<sub>3</sub>, and the 
corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For NO<sub>3</sub> radical, we employ polar versions of two dimensionless normal mode coordinates, Q<sub>3x</sub> and 
Q<sub>3y</sub> (degenerate components of in-plane asymmetric stretching vibration) to compute *ab initio* adiabatic PESs for five lowest 
electronic states and NACTs within them. The adiabatic potential energies are calculated in complete active space self-consistent field 
(CASSCF) level implementing 6-31G<sup>∗∗</sup> basis set and 9 electrons in 8 orbitals (9e, 8o) configuration active space (CAS). NACTs are 
computed using coupled-perturbed multi configuration space self-consistent field (CP-MCSCF) method implemented in MOLPRO quantum chemistry 
software. The first coordinate, &rho; ranges from 0.0 to 5.0 and the second coordinate, &phi; spans from 0 to 2&pi; such that total number of
grid points will be 50 X 180. Integration is performed along 'Path 6' and Numpy Binary format is used for both input and output. 

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_6'


>`adt num -mfile Tau_Rho.npy -nfile Tau_Phi.npy -efile Adiabatic_PES.npy -intpath 6 -nb`

---

Alternatively, user can employ the numerical segment of this package as a python module. Execution of python script file, 'script.py' generates
the output folder, 'ADT_script_2D_path6', which contains same results as 'ADT_numeric_6'.  
 
---
## Command:

>`python script.py`

---

