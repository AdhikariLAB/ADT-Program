import setuptools  # "magic" import
from numpy.distutils.core import setup, Extension
from setuptools import find_packages



# Build and Run or directly run with specified flags
# or run > python setup.py -h for details about the commands and flags

# By default this will install the package in root of the system and install an command line utility 'adtcli'
# if you dont have root privilage then install it with '--user' flag to install the package in a 
# local package folder for python(usuall ~/.local/lib/python2.7/site-packages/) and install the command line utility in
# respective folder (e.g /.local/bin/)
# Or you can just install this in a python virtualenv.




#uncomment one of the following as your requirement

#for default_fortran compiler(usually gfortran) without parallel
# python setup.py config_fc --fcompiler=intelem  install
# fort_args = []
# lib_links = []


# for ifort with openmp parallel flags
# python setup.py config_fc --fcompiler=intelem  install
# fort_args = ['-qopenmp']
# lib_links = ['-liomp5']



#for gfortran with openmp parallel flags
# python setup.py config_fc --fcompiler=gnu95  install
fort_args = ['-fopenmp']
lib_links = ['-lgomp']



lib = Extension(name='adtmod', 
            sources=['adt/numeric/nummod.f90']
            extra_f90_compile_args=fort_args,
            extra_link_args=lib_links)



setup(
    name='ADT_Program_Package',
    version='0.0.0',
    description='A description',
    author='whoever',
    author_email='whoever',
    classifiers=[
        'Development Status :: First Release',
        'Intended Audience :: Theoretical Chemistry',
        'Topic :: Software Development :: Build Tools',
        'License ::  MIT License',
        'Programming Language :: Python :: 2.7'
    ],
    zip_safe=True,
    install_requires=['numpy >= 1.10', 'h5py', 'six'],
    python_requires='>=2.7, <3',
    packages=find_packages(),
    ext_modules=[lib],
    entry_points={
        'console_scripts': [
            'adtcli = adt.adt:main',
        ],
    }
)
