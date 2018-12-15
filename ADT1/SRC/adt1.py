#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = 
ND=N*(N-1)/2

#formation of adiabatic potential energy matrix
U_Matrix = ModuleBase.adiabatic(N)

#writing of adiabatic potential energy matrix elements 
# fl = open('U_MATRIX.DAT','w')

# for i in range(1,N+1):
#   for j in range(1,N+1):
#     fl.write(' U_%r_%r = ' % (i,j))      
#     E = U_Matrix[i-1][j-1]  
#     fl.write(E)
#     fl.write('\n')
#     fl.write('\n')


txt =""
for i in range(1,N+1):
  for j in range(1,N+1):
    txt += ' U_%s_%s = ' % (i,j)
    txt += U_Matrix[i-1][j-1] +"\n\n"  

with open('U_MATRIX.DAT', "w") as f:
    f.write(txt)

info="""

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


#########################################################################

#The result file contains elements of Adiabatic Potential Energy Matrix.

#########################################################################


"""
with open('INFORMATION.DAT','w') as f:
    f.write(info)