#!/bin/bash

########################################################################################################################
#                                                                                                                      #
#    Execution of this script file is the initialization step of 'ADT' program package. This script compiles           #
#    the fortran file, adt.f90 to a python module file, adt_module.so using f2py. User can modify the flags of         #
#    f2py commands according to chosen compiler options. Moreover, openmp parallelization scheme is implemented        #
#    while solving the coupled adiabatic to diabatic transformation (ADT) equations. Hence, user has to specify        #
#    appropriate flags for parallel computation. In addition, the path of 'adt_final.py' is aliased as 'adt' in        #
#    the '.bashrc' file. Therefore, user can execute any symbolic manipulation or numerical calculation from any       #
#    directory using 'adt ... ... ...'  command. The associated arguments are thoroughly explained in 'adt_final.py'   #
#    file.                                                                                                             #
#                                                                                                                      #
#    To initialize 'ADT':                                                                                              #
#                                                                                                                      #
#         run ./init.sh in terminal (this step has to be performed only one time)                                      #
#                                                                                                                      #
#    Written by Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit     #
#    Adhikari                                                                                                          #   
#                                                                                                                      #
########################################################################################################################

#Creation of adt_module.so from adt.f90 using f2py command

# for gfortran with openmp parallelization
f2py -c adt.f90 -m adt_module --f90flags='-fopenmp' -lgomp only: get_angle amat 

########################################################################################################################
#                                                                                                                      #
# User can modify the above f2py command according to the compiler options and other flags. For detailed desciption of #
# 'f2py' command and its flags, visit 'https://docs.scipy.org/doc/numpy/f2py/'.                                        #
#                                                                                                                      #
# As for example, in case of ifort compiler, the 'f2py' command can be written as (for openmp parallelization):        #
# f2py -c adt.f90 -m adt_module --fcompiler=intelem --f90flags='-qopenmp' -liomp5 only: get_angle amat                 #
#                                                                                                                      #
########################################################################################################################

#Aliasing the path of 'adt_final.py' as adt

echo alias adt='"'python "'$(pwd)/adt_final.py'"  '"'>>~/.bashrc
source ~/.bashrc

########################################################################################################################
