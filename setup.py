__doc__='''
Script to install the 'ADT' package. 

Build and Run or directly run with specified flags.
execute python setup.py -h for details about the commands and flags

By default this will install the package in root of the system and install an command line utility 'adt'
if you dont have root privilage then install it with '--user' flag to install the package in a 
local package folder for python(usuall ~/.local/lib/python2.7/site-packages/) and install the command line utility in
respective folder (e.g /.local/bin/)
Or you can just install this in a python virtualenv.

'''

__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''
import setuptools
from numpy.distutils.core import setup, Extension
from setuptools import find_packages







#uncomment one of the following blocks as your requirement

#for default_fortran compiler(usually gfortran) without parallel
# python setup.py install
fort_args = []
lib_links = []


# for ifort with openmp parallel flags
# python setup.py config_fc --fcompiler=intelem  install
# fort_args = ['-qopenmp']
# lib_links = ['-liomp5']



#for gfortran with openmp parallel flags
# python setup.py config_fc --fcompiler=gnu95  install
#fort_args = ['-fopenmp']
#lib_links = ['-lgomp']



lib = Extension(name='adt.numeric.adtmod', 
            sources=['adt/numeric/nummod.f90'],
            extra_f90_compile_args=fort_args,
            extra_link_args=lib_links)



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
