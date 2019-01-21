#!/bin/bash


# with openmp parallelization, NOTE: 'fopenmp' is a gfortran flag
f2py -c adt.f90 -m adt_module --f90flags='-fopenmp' -lgomp only: get_angle amat 
# for ifort, following flags may not work universally
f2py -c adt.f90 -m adt_module --fcompiler=intelem --f90flags='-qopenmp' -liomp5 only: get_angle amat 
# f2py -c adt.f90 -m adt_module only: get_angle amat                                 #non-openmp
echo alias adt='"'python "'$(pwd)/adt_final.py'"  '"'>>~/.bashrc
source ~/.bashrc