#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase

N = 'NO OF STATE'
ND=N*(N-1)/2

#formation of adiabatic to diabatic transformation (ADT) matrix
A_Matrix = ModuleBase.matman(N)

#differentiating the relevant elements of ADT matrix
LHSelems = ModuleBase.elemgradselect(A_Matrix)

#formation of nonadiabatic coupling matrix (NACM)
NACM = ModuleBase.nacm(N)

#multiplication of negative of NACM and ADT matrix
R_Matrix =ModuleBase.multiply(ModuleBase.negative(NACM),A_Matrix)

#collecting the relevant elements of the above product matrix
RHSelems = ModuleBase.elemtauselect(R_Matrix)

#writing the partially substituted form of ADT equations
fl = open('ADT_EQUATIONS_PARTIAL.DAT','w')

count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    fl.write(' c(%r) = cos(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')
    fl.write(' s(%r) = sin(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')
    fl.write(' ic(%r) = sec(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')
    fl.write(' is(%r) = cosec(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')

fl.close()

Equations = ModuleBase.equation_partial(LHSelems,RHSelems,N)

#creating information file
fl1 = open('INFORMATION.DAT','w')
fl1.write('\n')
fl1.write('################################################################################################################################')
fl1.write('\n')
fl1.write('#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)')
fl1.write('\n')
fl1.write('################################################################################################################################')
fl1.write('\n')
fl1.write('\n')
fl1.write('######################################################################################################')
fl1.write('\n')
fl1.write('#The result file contains partially substituted Adiabatic to Diabatic Transformation (ADT) Equations.')
fl1.write('\n')
fl1.write('######################################################################################################')
fl1.write('\n')
fl1.write('\n')

fl1.close()
