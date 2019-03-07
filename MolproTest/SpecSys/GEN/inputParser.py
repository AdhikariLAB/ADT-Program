import numpy as np 


geomData = np.loadtxt('equiGeom.dat', 
    dtype={'names': ('names', 'mass', 'x','y','z'),'formats': ('S1', 'f8', 'f8','f8','f8')})
atomNames= geomData['names']
atomsMass= geomData['mass']
equiGeom = np.array([geomData[i] for i in ['x','y','z']])

freq = np.loadtxt('freq.dat')
wilsonMatrix = np.loadtxt('wilson.dat')



################################
#this section is to be continued using configparser
eInfo={
    'method':'multi',
    'basis':'6-31G**',
    'cas':'occ,13;closed,5',
    'wf':'wf,23,1,1,0',
    'state':'2'
}


nInfo={
    'method':'',
    'basis':'',
    'cas':'',
    'wf':'',
    'state':''
}
#####################################################################
state = 4

molproEquiTemplate='''
***, Molpro template created from ADT program for "WhatEver"
memory,100,m
file,2,molpro.wfu,new;

basis=6-31G**;

symmetry,nosym


geometry=geom.xyz

{{uhf;accu,5;orbital,2200.2}}
{{multi;occ,13;closed,5; wf,23,1,1,0;state,{state}; maxiter,40;orbital,2140.2;}}

show, energy
table, energy
save,molEquiEnr.res,new
{{table,____; noprint,heading}}
'''.format(state=state)




molproGridTemplate='''
***, Molpro template created from ADT program for "WhatEver"
memory,100,m
file,2,molpro.wfu,new;

basis = 6-31G**

geometry = geometry.xyz

!eInfo will be here
{{multi;occ,13;closed,5; wf,23,1,1,0;state,{state}; maxiter,40;start,2140.2;orbital,2140.2;}}

show, energy
table, energy
save,molenr.res,new
{{table,____; noprint,heading}}



'''.format(state=state)





nactPairs = [[i,j] for i in range(1,state+1) for j in range(i+1,state+1)]


nactTemp= "basis={}\n".format(eInfo['basis'])

for ind in range(0,len(nactPairs),5):
    nactTemp += "\n{{multi;occ,13;closed,5; wf,23,1,1,0;state,{};start,2140.2;\n".format(state)
    pairChunk = nactPairs[ind:ind+5]
    for i,j in pairChunk:
        nactTemp += "cpmcscf,nacm,{f}.1,{s}.1,accu=1.d-10,record=51{f}{s}.1;\n".format(f=i,s=j)
    nactTemp += "}\n"



for i,j in nactPairs:
    nactTemp +="""
force;nacm,51{f}{s}.1;varsav
table,gradx,grady,gradz;
save,ananac{f}{s}.res,new;
{{table,____;  noprint,heading}}

""".format(f=i,s=j)

molproGridTemplate +=nactTemp
print molproGridTemplate


