
# ADT PROGRAM
A generalised ADT program for analytical and numerical calculation. This is applicable for any
'N' electronic state sub-Hilbert space. The analytical segment can be used to generate symbolic
expressions of eigth adiabatic to diabatic transformation (ADT) quantities (elements of adiabatic
potential energy matrix, elements of nonadiabatic coupling matrix, ADT matrix elements, partially
substituted ADT equations, completely substituted ADT equations, elements of coefficient matrix of
gradient of ADT angles, elements  of coefficient matrix of nonadiabatic coupling terms (NACTs) and
the diabatic potential energy matrix elements) for any arbitrary number of coupled electronic
states. On the other hand, the numerical portion computes ADT angles, ADT matrix elements, diabatic
potential energy matrix elements and residue of ADT angles for any 'N' coupled electronic states
with multiple degrees of freedom. Any user can solve the differential equations along eight
different paths over the nuclear configuration space (CS).


## *Authors:*

Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

---

## Requirements: 
1. Fortran compiler
2. Python 2.7  

### Python dependencies required:
1. NumPy >= 1.13.0
2. h5py (only if you need HDF5 file I/O)

---

## Installation: 
For installing the package with the default parameters just run  

`python setup.py install`  

This will install the package in the root of the system, so you may need to run this with admin privilages. If you don't have root privilages or don't want to install the package in the home of the local user site-packages folder for python by running  

`python setup.py install --user`  

On successfull installation of the package, a command line utility `adt` will be installed which will be used to run the package directly from terminal.
A quick help about the installation can be found by running 
`python setup.py -h`
To know more about installation and building distribution please refer to the user manual of the package.   
### Installing with parallelization support:
By default the package will be installed using the default fortran compiler and without any parallelization support. OpenMP parallelization is implemented in this package and to take benifit of it user has to pass proper compiler flag during the installation.

To install with parallelization uing gfortran user has to run  
`python setup.py config_fc --fcompiler=gnu95  install`  
while to use Intel fortran compiler run  
`python setup.py config_fc --fcompiler=intelem  install`  
In addition to this, to properly link the OpenMP libraries user has to set the `fort_args` and `lib_links` variable, appropriate for the compiler, inside the setup.py script. Those variable, for gfortran and ifort are availabe in the script, and just uncomment the relavant lines to set them. For any other fortran compiler, user themselves has to find correct linker flags appropriate for that compiler.


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
│   ├── numeric/
│   │   ├── __init__.py
│   │   ├── adt_numeric.py          # script for numerical calculation
│   │   └── nummod.f90              # collection of necessary fortran subroutines
│   └── molpro/
│       ├── __init__.py
│       └── adt_molpro.py           # script for ab initio calculations by MOLPRO
├── test_runs/                      # folder containing sample calculations
├── license                         # icense information
├── user_manual.pdf                 # detailed instructions for users
└── readme.md                       # readme file
```



## Usage:

All the operation of the package can be done using the command line utility `adt`. Detailed usage of the command can be seen by running 
`adt -h`

The main command `adt` has three subcommands for three different type of calculation namely  
1. __ana__ : Calculate different [analytic expressions](#analytic) for by passing the subcommand `ana`.
2. __num__ : Subcommand `num` can be used to calculate different [numerical quantities](numerical).
3. __mol__ : Use subcommand 'mol' to calculate _ab-initio_ PESs and NACTs for a molecular species using [MOLPRO](#molpro) and subsequently calculate the numerical quantities.

Description for each of the segments are given in the following and can be be found in detail in the user manual. At any step of using this package user can see the help menu by passing the help flag `-h` or `--help`. During the runtime of any of the section the, all the relavant information and progress is saved in a logfile named 'ADT.log'.

### Analytic:
This section is employed to derive the following eight analytical quantities for any number of electronic states constituting the sub-Hilbert space,

1. Adiabatic potential energy matrix.
2. Non-Adiabatic Coupling Matrix (NACT).
3. ADT Matrix.
4. Partilly substituted ADT Equations.
5. Complete form of ADT equations.
5. Coefficient matrix of gradient of ADT angle.
6. Coefficient matrix of nact.
7. Diabatic potential energy matrix

The Analytical section is accessed using the subcommand `ana`. The manual for this section can be seen by passing the help flag `-h`,

```

$ adt ana -h
usage: adt ana [-h] -nstate NSTATE [-anajob [1-8]]

Devise analytical expressions of any one of the ADT quantities for a given number of states

optional arguments:
  -h, --help      show this help message and exit
  -anajob [1-8]   Specify the type of expression (default: 5 - completely substituted ADT equation) 

Required arguments:
  -nstate NSTATE  Number of states
```
The required argument `-nstate` takes the number of electronic states to derive the analytical expression for. On the other hand `-anajob` takes a numer ranging from 1 to 8 representing the type of type of expression listed above.


### Numerical: 
This segment takes the adiabatic potential energy values and non adiabatic coupling matrix elements and  calculates the ADT angles, ADT matrices, diabatic potential energy matrices and ADT angle residues at any number of grid points and for N-dimensional sub-Hilbert space 
(N is arbitrary) in the interested domain of nuclear configuration space (CS).

Following is the help menu for the numerical section,
```
$ adt num -h
usage: adt num [-h] -nfile1 FILE -nfile2 FILE [-intpath [1-8]] [-efile FILE]
               [-nstate NSTATE] [-ofile FILE] [-n N] [-h5] [-nb] [-txt]

Calculate numerical results for a given number of electronic states along a specific path

optional arguments:
  -h, --help      show this help message and exit
  -intpath [1-8]  Specify the path for calculation (default: 1).
                   
  -efile FILE     Specify the Energy file for calculating the diabatic potential energy matrix elements.
                   
  -nstate NSTATE  Specify the number of states to do the calculation.
                  By default it includes all the data for calculation.
                    
  -ofile FILE     Specify the output file name (w/o extension) (default: 'ADT_numeric').
                   
  -n N            Specify number of OpenMP threads to use for parallel calculation. 
                  Applicable only when installed using OpenMP support.
                  (default: Maximum avilable threads)
                   
  -h5             Write results in a HDF5 file (.h5). 
                  Fast IO, smaller file size and hierarchical filesystem-like data format,
                  preferable for saving and sharing large datasets in an organised way.
                   
  -nb             Write results in Numpy binary file (.npy). 
                  Preferable when working with numpy for its much faster IO and easy portability.
                   
  -txt            Write results in a text file. (default behaviour).

Required arguments:
  -nfile1 FILE    Specify the NACT file along first coordinate.
                   
  -nfile2 FILE    Specify the NACT file along second coordinate.
```

The two required arguments `-nfile` and `-nfile2` takes the NACT filenames along the first and second coordinate respectively while the energy file provided using `-efile` is only needed when user need to calculate the diabatic potential energy matrix also. The output files are saved in a folder for Numpy binary format and Plain txt while for HDF5 all the output data is saved in a single file. Name of the file/folder is set using `-ofile` argument. Any of three file types i.e. HDF5, Numpy binary format and Plain txt and be used as file I/O for this package and can be set through `-h5`,`-nb` and `-txt` respectively. Though by default calculations are done taking all the sates provided in the input data, user can use the argument `-nstate` to set a differnet number of states for calculation. When the package is installed using OpenMP support, the program will run with maximum OpwnMP threads availble, or the user can set the number of threads using the argument `-n`.

__Some key things to note here:__ 
* The coordinates of grid points can be chosen from any coordinate system like 
    cartesian coordinate, plane polar coorinate, hyperspherical coordinate etc.
* ADT matrix is formed by multiplying the elementary rotation matrices in a 
    particular order. they are arranged in increasing order of their second indeces
    and then in ascending order of their first indeces keeping the second index fixed. 
    [for four electronic states, the order is a(1,2)*a(1,3)*a(2,3)*a(1,4)*a(2,4)*a(3,4)] 
* The ADT equations are numerically solved by 8th order Runge-Kutta method.
* Numerical integration of the adt equations can be solved along infinite number of
    paths in a two-dimensional coordinate system. here the differential equations can 
    be solved along eight (8) possible paths taking the order and direction of the two coordinate 

* The magnitude of adt angle residue will be meaningful only if the second coordinate
     (supplied by the user) is bound that is it must form a closed contour. 







### Molpro:
In addition to providing the energy and nonadiabatic coupling terms through the input files, user can directly calculate the same using this section of the package providing the information about the molecular species.

Help message can also be seen for this section by passing the help flag `-h`.
```
$ adt mol -h
usage: adt mol [-h] -sys SYS [-config FILE] [-atomfile FILE] [-geomfile FILE]
               [-freqfile FILE] [-wilsonfile FILE] [-intpath [1-8]]
               [-ofile FILE] [-n N] [-h5] [-nb] [-txt]

Calculate the Adiabatic potential energy surfaces(PESs) and noadiabatic coupling matrix(NACM)
and subsequently calculate the Numerical quantities using the numerical section.

optional arguments:
  -h, --help        show this help message and exit
  -config FILE      Specify the molpro configuration file containing
                    the necessary keywords of MOLPRO software.
                     (default: molpro.config). 
                     
  -atomfile FILE    Specify the information file constituting atomic
                    symbols and atomic masses
                    (default: atomfile.dat). 
                     
  -geomfile FILE    Specify the geometry file containing the initial
                    grid point in "xyz" format (in Angstrom)
                    (default: geomfile.dat). 
                    (Igonred for scattering system). 
                     
  -freqfile FILE    Specify the frequency information file, where
                    frequencies of normal modes are written in cm-1
                    (default: frequency.dat). 
                    (Igonred for scattering system). 
                     
  -wilsonfile FILE  Specify the filename containing the Wilson matrix
                    of a molecular species (default: wilson.dat).
                    (default: wilson.dat).
                    (Igonred for scattering system). 
                     
  -intpath [1-8]    Specify the path for calculation (default: 1).
                     
  -ofile FILE       Specify the output file name (w/o extension) (default: 'ADT_numeric').
                     
  -n N              Specify number of OpenMP threads to use for parallel calculation. 
                    Applicable only when installed using OpenMP support.
                    (default: Maximum avilable threads)
                     
  -h5               Write results in a HDF5 file (.h5).
                    Fast IO, smaller file size and hierarchical filesystem-like data format,
                    preferable for saving and sharing large datasets in an organised way.
                     
  -nb               Write results in Numpy binary file (.npy). 
                    Preferable when working with numpy for its much faster IO and easy portability.
                     
  -txt              Write results in a text file. (default behaviour).

Required arguments:
  -sys SYS          Specify type of the molecular process 
                    (molecular species for spectroscopic calculation or reactive moeities for scattering calculation) 
                    (Available options: spectroscopic or scattering)

```

Most of the options namely `intpath`,`-ofile`,`-n`,`-nb`,`-txt` are exactly like the previous numerical section. Using `-sys` argument user has to define the type of the system you want to do _ab-initio_ of i.e. spectroscopic or scattering, while `-config` is the argument where user has to provide  filename containing the configuration to run the MOLPRO quantum chemistry package. Details about the system, under study, is provided using the arguments `-atomfile`,`-geomfile`,`-freqfile` and `-wilsonfile` which takes the information about atoms, equilibrium geometry, frequency of the normal modes and wilson matrix respectively, where the last three argument is only needed for a spectroscopic system. Details about the structure can be found in the user manual.


## Test Runs:
Some examples  for each of the three section are included in the 'Test_Runs' folder to easily understand the workflow of the package. 
### Files and folder structure for the example 'Test_Runs' folder:
```

```