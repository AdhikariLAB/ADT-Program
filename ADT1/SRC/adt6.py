#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase

N = 'NO OF STATE'
ND=N*(N-1)/2

#formation of adiabatic to diabatic transformation (ADT) matrix
A_Matrix = ModuleBase.matman(N)

#collecting the relevant elements of ADT matrix
Aelems = ModuleBase.elemselect(A_Matrix)

#differentiating the collected elements of ADT matrix
LHSelems = []

for i in range(len(Aelems)):
    LHSelems.append(ModuleBase.diff(Aelems[i],N))

#writing the coefficient matrix of the gradient of the ADT angles
fl = open('GRADCOEFF.DAT','w')

count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    fl.write(' c(%r) = cos(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')
    fl.write(' s(%r) = sin(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')

fl.close()

gradrow = 0 
for i in LHSelems:
    gradrow += 1
    ModuleBase.extractgrad(i,N,gradrow)

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
fl1.write('##########################################################################################################################')
fl1.write('##############')
fl1.write('\n')
fl1.write('#The result file contains elements of Coefficient Matrix corresponding to Gradient of Adiabatic to Diabatic Transformation') 
fl1.write('(ADT) Angles.')
fl1.write('\n')
fl1.write('##########################################################################################################################')
fl1.write('##############')
fl1.write('\n')
fl1.write('\n')

fl1.close()
