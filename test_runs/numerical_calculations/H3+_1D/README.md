

## Description


This folder contains adiabatic potential energy surfaces (PESs) and nonadiabatic coupling terms (NACTs) for H<sub>3</sub><sup>+</sup>, and the
corresponding ADT quantities (ADT angles, residue of ADT angles, ADT matrices and diabatic potential energy matrices) calculated by this 
program package. For 1D ADT calculation of H<sub>3</sub><sup>+</sup>, we employ Jacobi coordinate system, where &r;, &R;, &gamma; and &q; 
are chosen as 2 Bohr, 1.732050808 Bohr, 90<sup>o</sup> and 0.5 Bohr, respectively. We use other coordinate, &phi; to compute *ab initio* 
adiabatic PESs for three lowest electronic states and NACTs within them. The adiabatic potential energies are calculated in multireference 
configuration interaction (MRCI) level implementing Dunning’s correlation consistent quintuple zeta (cc-pV5Z) basis set and 2 electrons in 
10 orbitals (2e, 10o) configuration active space (CAS). NACTs are computed using coupled-perturbed multi-configuration space self-consistent 
field (CP-MCSCF) method implemented in MOLPRO quantum chemistry software. The circular coordinate, &phi; spans from 0 to 2&pi; such that 
total number of grid points wil be 121. Integration is performed along this circular coordinate and simple text format is used for both input and output.

---
## Command:

The following command is used to generate the results, which are stored in the folder 'ADT_numeric_1D'


>`adt num -nfile Tau_Phi.dat -efile Adiabatic_PES.dat`

---

Alternatively, user can employ the numerical segment of this package as a python module. Execution of python script file, 'script.py' generates
the output folder, 'ADT_script_1D', which contains same results as 'ADT_numeric_1D'.  
 
---
## Command:

>`python script.py`

---

