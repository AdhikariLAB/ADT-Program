#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = int(sys.argv[1])
ND=N*(N-1)/2

#formation of adiabatic to diabatic transformation (ADT) matrix
A_Matrix = ModuleBase.matman(N)

#formation of nonadiabatic coupling matrix (NACM)
NACM = ModuleBase.nacm(N)

#multiplication of negative of NACM and ADT matrix
R_Matrix = ModuleBase.multiply(ModuleBase.negative(NACM),A_Matrix)

#collecting the relevant elements of the above product matrix
RHSelems = ModuleBase.elemselect(R_Matrix)

#writing the coefficient matrix of the nonadiabatic coupling terms (NACT)
txt= ""
count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
        .format(c=count, i1=index1, i2=index2)

with open('TAUCOEFF.DAT','w') as f:
    f.write(txt)

# taurow = 0 
# for i in RHSelems:
#     taurow += 1
#     ModuleBase.extracttau(i,N,taurow)

for j,i in enumerate(RHSelems,1):
    ModuleBase.extracttau(i,N,j)

#creating information file
info="""

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


###############################################################################################################

#The result file contains elements of Coefficient Matrix corresponding to Nonadiabatic Coupling Terms (NACTs).

###############################################################################################################

"""
with open('INFORMATION.DAT','w') as f:
    f.write(info)
