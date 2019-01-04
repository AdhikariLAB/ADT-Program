#!/bin/bash

f2py -c adt.f90 -m adt_module only: get_angle amat
echo alias adt='"'python "'$(pwd)/adt_final.py'"  '"'>>~/.bashrc
source ~/.bashrc