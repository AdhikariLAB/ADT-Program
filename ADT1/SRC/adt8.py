#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase

N = 'NO OF STATE'
ND=N*(N-1)/2

#formation of adiabatic potential energy matrix
A_Matrix = ModuleBase.matman(N)

#formation of diabatic potential energy matrix
W_Matrix = ModuleBase.diabatic(A_Matrix)

#writing of diabatic potential energy matrix elements
fl = open('W_MATRIX.DAT','w')

count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    fl.write(' c(%r) = cos(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')
    fl.write(' s(%r) = sin(theta(%r,%r))\n' % (count,index2,index1))
    fl.write('\n')

for i in range(1,N+1):
    for j in range(1,N+1):
        fl.write(' W_%r_%r = \n' % (i,j))      
        E = W_Matrix[i-1][j-1]  
        length = len(E)
        line = length/140                                                                                                                          
        for k1 in range(1,line+1):
          start = 0 + (k1-1)*140
          end = k1*140 
          segment = E[start:end]
          fl.write(' ')
          fl.write(segment)
          fl.write('\n')
        fl.write(' ')
        fl.write(E[line*140:])
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
fl1.write('# The result file contains elements of Diabatic Potential Energy Matrix.')
fl1.write('\n')
fl1.write('#########################################################################')
fl1.write('\n')
fl1.write('\n')

fl1.close()
