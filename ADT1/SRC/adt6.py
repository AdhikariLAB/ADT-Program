#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = int(sys.argv[1])
ND=N*(N-1)/2

#formation of adiabatic to diabatic transformation (ADT) matrix
A_Matrix = ModuleBase.matman(N)

#collecting the relevant elements of ADT matrix
Aelems = ModuleBase.elemselect(A_Matrix)

#differentiating the collected elements of ADT matrix
# LHSelems = []
# for i in range(len(Aelems)):
#     LHSelems.append(ModuleBase.diff(Aelems[i],N))
LHSelems = [ModuleBase.diff(el,N) for el in Aelems]

#writing the coefficient matrix of the gradient of the ADT angles

txt= ""
count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
        .format(c=count, i1=index1, i2=index2)

with open('GRADCOEFF.DAT','w') as f:
    f.write(txt)



# gradrow = 0 
# for i in LHSelems:
#     gradrow += 1
#     ModuleBase.extractgrad(i,N,gradrow)

for j,i in enumerate(LHSelems,1):
    ModuleBase.extractgrad(i,N,j)
#creating information file

info = """

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


########################################################################################################################################

#The result file contains elements of Coefficient Matrix corresponding to Gradient of Adiabatic to Diabatic Transformation
(ADT) Angles.

########################################################################################################################################

"""


with open('INFORMATION.DAT','w') as f:
    f.write(info)