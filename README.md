
# ADT PROGRAM
A generalised ADT program for analytical and numerical calculation. This is applicable for any
'N' electronic state sub-Hilbert space. The analytical segment can be used to generate symbolic
expressions of six quantities related to adiabatic to diabatic transformation (ADT) and those 
are: ADT matrix elements, partially substituted ADT equations, completely substituted ADT equations, 
elements of coefficient matrix of gradient of ADT angles, elements  of coefficient matrix of 
nonadiabatic coupling terms (NACTs) and the diabatic potential energy matrix elements for any 
arbitrary number of coupled electronic states. On the other hand, the numerical portion computes 
ADT angles, ADT matrix elements, diabatic potential energy matrix elements and residue of ADT 
angles for any 'N' coupled electronic states with multiple degrees of freedom. An user can solve 
the differential equations along any one of eight different paths over the nuclear configuration 
space (CS).


## *Authors:*

Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Satyam Ravi, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari

---

## Requirements: 
1. Fortran compiler
2. Python 2.7 or 3.x  

### Python dependencies required:
1. setuptools
2. NumPy >= 1.13.0
3. h5py (only if you need HDF5 file I/O)

---

## Installation: 

#### Setting Fortran compiler for installation

The 'ADT' program uses OpenMP directives in the Fortran code and to properly install the package, the code has to be compiled correctly with appropriate flags. To do so one has to let the build system know what Fortran compiler to use, during installation. The simplest way to do this is to set the `F90` environment variable before installation. 
To install using `gfortran` compiler run (in Bash)

```bash
export F90=gfortran
```
or to use Intel Fortran compiler:

```bash
export F90=ifort
```
If no `F90` environment variable is specified then the program is installed using `gfortran` compiler by default. To use other compilers or to fine tune the installation process, modify the `setup.py` script.

#### Installing 'ADT'
For installing the package with the default parameters just run:  
```
python setup.py install
```

This will install the package, usually in the root of the system, so you may need to run this with admin privilages. If you don't have root privilages or don't want to install the package in the root of the system, user can install this package in home of the local user site-packages folder for python by 
running  

```
python setup.py install --user
```

Alternatively, if `pip` is available, user can install the package by executing the following command inside the source directory

```
pip install .
```

On successfull installation of the package, a command line utility `adt` will be created, which will be used to run the package directly from terminal.

A quick help about the installation can be found by running `python setup.py -h`  

To know more about installation and building distribution, please refer to the user manual of the package.   


---

## Files and Folder Structure of the package:


```bash
ADT/
├── setup.py                        # Package Installer script
├── adt/
│   ├── adt.py                      # central script to run the package
│   ├── __init__.py
│   ├── analytic/
│   │   ├── __init__.py
│   │   ├── adt_analytic.py         # script for analytic formulation
│   │   └── anamod.py               # collection of necessary python functions
│   │
│   ├── numeric/
│   │   ├── __init__.py
│   │   ├── adt_numeric.py          # script for numerical calculation
│   │   └── nummod.f90              # collection of necessary fortran subroutines
│   │
│   ├── optimization/
│   │   ├── __init__.py
│   │   └── optimize.py             # script for geometry optimization, frequency and wilson matrix by MOLPRO or Gaussian or 
│   │                               # Gamess package   
│   └── molpro/
│       ├── __init__.py
│       └── adt_molpro.py           # script for ab initio PESs and NACTs calculations by MOLPRO  
│ 
├── test_runs/                      # folder containing sample calculations
├── LICENSE                         # license information
├── user_manual.pdf                 # detailed instructions for users
└── README.md                       # readme file
```

## Usage:

All the operation of the package can be done using the command line utility `adt`. The main command `adt` has three subcommands for three different 
types of calculation, namely,
  
1. __ana__ : Calculate different [analytic expressions](#analytic) by using the subcommand `ana`.
2. __num__ : Subcommand `num` can be used to calculate different [numerical quantities](#numerical).
3. __opt__ : The `opt` subcommand is used to calculate [optimized geometry](#optimization), frequencies and wilson matrix of a Spectroscopic system, required for the _ab-initio_.
4. __mol__ : Use subcommand 'mol' to calculate _ab-initio_ PESs and NACTs for a molecular species using [MOLPRO](#molpro) and subsequently, calculate the numerical quantities

Description for each of the above segments can be found in detail in the user manual. At any step of using this package, user can see the help menu by using the help flag `-h` or `--help`. During the runtime of any of the above sections, all the relevant information and progress is saved in a logfile named 'ADT.log'.


### Analytic:
This section is employed to derive the following six analytical quantities for any number of electronic states constituting the sub-Hilbert space:

1. ADT Matrix.
2. Partilly substituted ADT Equations.
3. Complete form of ADT equations.
4. Coefficient matrix of gradient of ADT angle.
5. Coefficient matrix of NACT.
6. Diabatic potential energy matrix

The Analytical section is accessed using the subcommand `ana`. The manual for this section can be seen by using the help flag `-h`,


### Numerical: 
This segment takes the values of adiabatic potential energy surfaces and non adiabatic coupling matrix elements of N-dimensional sub-Hilbert 
space for a set of nuclear grid points, and calculates the ADT angles, ADT matrices, diabatic potential energy matrices and ADT angle 
residues for the same set of grid points.

__Some key things to note here:__ 

* The coordinates of grid points can be chosen from any coordinate system like 
    normal mode coordinate, hyperspherical coordinate and Jacobi coordinate.

* ADT matrix is formed by multiplying the elementary rotation matrices in any order. 

* The ADT equations are numerically solved by 8th order Runge-Kutta method.

* Numerical integration of the adt equations can be solved along infinite number of
    paths in a two-dimensional coordinate system. here the differential equations can 
    be solved along eight (8) possible paths taking the order and direction of the two coordinate 

* The magnitude of adt angle residue will be meaningful only if the second coordinate
     (supplied by the user) is bound that is it must form a closed contour. 

### Optimization:
Perform geometry optimization to calculate optimized geometry, frequencies and wilson matrix of a Spectroscopic system, required for the _ab-initio_. Presently Molpro Gamess and Gamess can be used for this purpose.

### Molpro:
If the adiabatic PESs and NACTs are not available, user can directly calculate those interfacing the MOLPRO providing the information about the molecular species

## Citation:
This work is published in [Journal of Chemical Theory and Computation](https://pubs.acs.org/journal/jctcce) in paper :  
**ADT : A Generalized Algorithm and Program for Beyond Born-Oppenheimer Equations of 'N' Dimensional Sub-Hilbert Space** (DOI : [https://doi.org/10.1021/acs.jctc.9b00948](https://doi.org/10.1021/acs.jctc.9b00948))  
Consider citing this paper if you use this package 
