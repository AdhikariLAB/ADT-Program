#!/bin/bash

f2py -c adt.f90 -m adt_module --f90flags='-fopenmp' -lgomp only: get_angle amat  #with openmp parallelization
# f2py -c adt.f90 -m adt_module only: get_angle amat                                 #non-openmp
echo alias adt='"'python "'$(pwd)/adt_final.py'"  '"'>>~/.bashrc
source ~/.bashrc