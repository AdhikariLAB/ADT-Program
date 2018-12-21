
import sys
import subprocess

from sys import argv
script, path=argv

fl = open('input.dat','r')

param = fl.readline().strip()

paramlist = param.split(',')

paramdict = {}

for i in range(len(paramlist)):
  paramdict[paramlist[i].split('=')[0].strip()]=paramlist[i].split('=')[1].strip()

fl.close()

print paramdict

print script

if [paramdict['jobtype'] == 'analytic']:
   job = 1
   subprocess.call(["touch","%sadt.f" % path])
   print job 
