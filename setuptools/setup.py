import setuptools  # "magic" import
from numpy.distutils.core import setup, Extension


lib = Extension(name='adt.numeric.adtmod', sources=[
                'Working/numeric/nummod.f90'])

setup(
    name='ADT_Program_Package',
    version='0.0.0',
    description='A description',
    author='whoever',
    author_email='whoever',
    classifiers=[
        'Development Status :: First Release',
        'Intended Audience :: Theoretical Chemistry People',
        'Topic :: Software Development :: Build Tools',
        'License ::  MIT License',
        'Programming Language :: Python :: 2.7'
    ],
    install_requires=['numpy >= 1.10', 'h5py'],
    python_requires='>=2.7, <3',
    packages=['adt'],
    package_dir={'adt': 'Working'},
    ext_modules=[lib],
    entry_points={
        'console_scripts': [
            'adtcli = adt.adt:main',
        ],
    }
)
