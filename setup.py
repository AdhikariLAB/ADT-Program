__doc__='''
Script to install the 'ADT' package. 

Build and Run or directly run with specified flags.
execute python setup.py -h for details about the commands and flags

By default this will install the package in root of the system and install an command line utility 'adt'
if you dont have root privilage then install it with '--user' flag to install the package in a 
local package folder for python(usuall ~/.local/lib/pythonX.X/site-packages/) and install the command line utility in
respective folder (e.g /.local/bin/)

'''

__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Satyam Ravi, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''
import os
import setuptools
from numpy.distutils.core import setup, Extension
from setuptools import find_packages



# During the installation the fortran code is comipled and as OpenMP is used here, proper compilation
# flags has to be used to compile the code for a particular fortran compiler. 

# By default the program uses gfortran as its fortran compiler. One can set the fortran compiler to be used through the environment varibale `F90`
# To use `gfortran` run in bash `export F90=gfortran` (or just clear the variable as `gfortran` is default)
# To use `ifort` run in bash `export F90=ifort`

# The program can only handle, automatically, `gfortran`,`ifort` and PGI fortran compilers (see below), for any other compiler user has to manully
# configure this script by setting `f90_flags` and `omp_lib` for that particular compiler and running the script by
# `python setup.py config_fc --fcompiler=<compiler_name>  install`
# for compiler_name, one can see the available  fortran compiler for the platform by running
# `python setup.py config_fc --help-fcompiler`



F90 = os.getenv("F90")

# if no envirnment variable sent then use gfortran by default
if F90 == None or F90 == "":
    F90 = 'gfortran'
    print("NOTE: No 'F90' environment variable set. Using 'gfortran' by default")

if F90 == "ifort":                      # for `ifort`
    f90_flags = ["-qopenmp","-O3"]
    omp_lib = ["-liomp5"]

elif F90 == "gfortran":                 # for `gfortran`
    f90_flags = ["-fopenmp", "-fPIC", "-O3"]
    omp_lib = ["-lgomp"]

elif F90 in ["pgfortran", "pgf90", "pgf95"]: # PGI fortran compilers
    f90_flags = ["-mp"]
    omp_lib = [""]

else:
    ll = "Environment variable 'F90={}' not recognized.\n Configure 'setup.py' manually to set up proper OpenMP library links".format(F90)
    raise RuntimeError( ll )




lib = Extension(name='adt.numeric.adtmod', 
            sources=['adt/numeric/nummod.f90'],
            extra_f90_compile_args=f90_flags,
            extra_link_args=omp_lib)



setup(
    name='adt',
    version='0.1.0',
    description='Installer for the ADT_Program_Package.',
    author='Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari',
    author_email='pcsa@iacs.res.in',
    license = "GNU GPLv2",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: CPython'
    ],
    keywords='Quantum chemistry, PES, NACT',
    project_urls={'Source Code':'https://github.com/AdhikariLAB/ADT-Program'},
    zip_safe=True,
    setup_requires=['numpy >=1.13.0'],
    install_requires=['numpy >=1.13.0'],
    extras_require={
        'h5':  ["h5py"]
    },
    python_requires='>=2.7',
    packages=find_packages(),
    ext_modules=[lib],
    entry_points={
        'console_scripts': [
            'adt = adt.adt:main',
        ],
    }
)
