#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase

N = 'NO OF STATE'
ND=N*(N-1)/2

#formation of adiabatic potential energy matrix
U_Matrix = ModuleBase.adiabatic(N)

#writing of adiabatic potential energy matrix elements 
fl = open('U_MATRIX.DAT','w')

for i in range(1,N+1):
  for j in range(1,N+1):
    fl.write(' U_%r_%r = ' % (i,j))      
    E = U_Matrix[i-1][j-1]  
    fl.write(E)
    fl.write('\n')
    fl.write('\n')

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
fl1.write('#########################################################################')
fl1.write('\n')
fl1.write('#The result file contains elements of Adiabatic Potential Energy Matrix.')
fl1.write('\n')
fl1.write('#########################################################################')
fl1.write('\n')
fl1.write('\n')

fl1.close()
