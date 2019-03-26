import setuptools  # "magic" import
from numpy.distutils.core import setup, Extension


lib = Extension(name='adtp.adtmod', sources=['Working/adt.f90'])

setup(
    name = 'adtpackage2',
    packages = ['adtp'],
    package_dir = {'adtp':'Working'},
    ext_modules = [lib],
    entry_points={
        'console_scripts': [
            'adtcli = adtp.adtmain:main',
        ],
    }
)