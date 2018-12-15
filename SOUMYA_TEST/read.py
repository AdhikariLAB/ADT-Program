#!/bin/python

import os

fl = open('input.dat','r')

param = fl.readline().strip()

paramlist = param.split(',')

paramdict = {}

for i in range(len(paramlist)):
  paramdict[paramlist[i].split('=')[0].strip()]=paramlist[i].split('=')[1].strip()

fl.close()

print paramdict

if [paramdict['jobtype'] == 'analytic']:
   job = 1
   os.cd('ADT$job')
   print job 
