#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

import ModuleBase, sys

N = int(sys.argv[1])
ND=N*(N-1)/2

#formation of nonadiabatic coupling matrix (NACM)
NACM = ModuleBase.nacm(N)
#writing of NAC matrix elements
# fl = open('NACM.DAT','w')

# for i in range(1,N+1):
#   for j in range(1,N+1):
#     fl.write(' NACM_%r_%r = ' % (i,j)) 
#     fl.write(NACM[i-1][j-1])     
#     fl.write('\n')
#     fl.write('\n')


txt =""
for i in range(1,N+1):
  for j in range(1,N+1):
    txt += ' NACM_%s_%s = ' % (i,j) 
    txt += NACM[i-1][j-1] + "\n\n"  

with open('NACM.DAT','w') as f:
    f.write(txt)


#creating information file
info= """

################################################################################################################################

#Authors are Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)

################################################################################################################################


###########################################################################

#The result file contains elements of Nonadiabatic Coupling Matrix (NACM).

###########################################################################


"""

with open('INFORMATION.DAT','w') as f:
    f.write(info)