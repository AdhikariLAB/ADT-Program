#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = int(sys.argv[1])
ND=N*(N-1)/2

#formation of adiabatic to diabatic transfromation (ADT) matrix 
A_Matrix = ModuleBase.matman(N)


txt= ""
count = 0
for index1 in range(2,N+1):
  for index2 in range(1,index1):
    count += 1
    txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
        .format(c=count, i1=index1, i2=index2)


for i in range(1,N+1):
  for j in range(1,N+1):
    txt += ' A_%s_%s = \n' % (i,j)   
    E = A_Matrix[i-1][j-1]  
    txt += "\n".join(E[i:i+140] for i in range(0,len(E),140))+"\n\n"

with open('A_MATRIX.DAT','w') as f:
    f.write(txt)


#creating information file
info ="""

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


###################################################################################

#The result file contains elements of Adiabatic to Diabatic Transformation Matrix.

###################################################################################


"""
with open('INFORMATION.DAT','w') as f:
    f.write(info)