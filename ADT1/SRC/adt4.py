#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = int(sys.argv[1])
ND=N*(N-1)/2

#formation of adiabatic to diabatic transformation (ADT) matrix
A_Matrix = ModuleBase.matman(N)

#differentiating the relevant elements of ADT matrix
LHSelems = ModuleBase.elemgradselect(A_Matrix)

#formation of nonadiabatic coupling matrix (NACM)
NACM = ModuleBase.nacm(N)

#multiplication of negative of NACM and ADT matrix
R_Matrix = ModuleBase.multiply(ModuleBase.negative(NACM),A_Matrix)

#collecting the relevant elements of the above product matrix
RHSelems = ModuleBase.elemtauselect(R_Matrix)

#writing the partially substituted form of ADT equations




txt = ""
count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    txt +="""
c({c}) = cos(theta({i1},{i2}))

s({c}) = sin(theta({i1},{i2}))

ic({c})= sec(theta({i1},{i2}))

is({c})= cosec(theta({i1},{i2}))
""".format(c=count, i1=index1, i2=index2)


with open("ADT_EQUATIONS_PARTIAL.DAT", "w") as f:
    f.write(txt)


Equations = ModuleBase.equation_partial(LHSelems,RHSelems,N)

#creating information file
info ="""

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


###################################################################################

#The result file contains partially substituted Adiabatic to Diabatic Transformation (ADT) Equations.

###################################################################################


"""
with open('INFORMATION.DAT','w') as f:
    f.write(info)