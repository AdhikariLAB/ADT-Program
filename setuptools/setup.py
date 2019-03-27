import setuptools  # "magic" import
from numpy.distutils.core import setup, Extension
from setuptools import find_packages

lib = Extension(name='adtmod', 
            sources=['adt/numeric/nummod.f90'],
            extra_f90_compile_args=['-fopenmp','-lgomp'],
            libraries=['gomp'])

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
