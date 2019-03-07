import numpy as np
import ConfigParser ,textwrap, sys, os


angtobohr = 1.88973
hbar = 0.063508
cInvToTauInv = 0.001883651


def parseConfig(conFigFile='molpro.config'):
    scf = ConfigParser.SafeConfigParser()
    scf.read(conFigFile)
    eInfo = dict(scf.items('eInfo'))
    nInfo = dict(scf.items('nInfo'))
    gInfo = dict(scf.items('gInfo'))
    mInfo = dict(scf.items('mInfo'))
    vModes = [int(i)-1 for i in mInfo['varying'].split(',')]
    return [eInfo, nInfo,vModes, gInfo]



def createTemplate(eInfo, nInfo):
    molproEquiTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        memory,100,m
        file,2,molpro.wfu,new;

        basis=6-31G**;

        symmetry,nosym


        geometry=geom.xyz

        {{uhf;accu,5;orbital,2200.2}}
        {{{method};{cas}; wf,{wf};state,{state};maxiter,40; orbital,2140.2;}}

        show, energy
        table, energy
        save,equienr.res,new
        {{table,____; noprint,heading}}

        ---

        '''.format(method=eInfo['method'],
                    state=eInfo['state'],
                    wf =  eInfo['wf'],
                    cas = eInfo['cas']))

    with open('equi.com' ,'w') as f:
        f.write(molproEquiTemplate)



    molproGridTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        memory,100,m
        file,2,molpro.wfu,new;

        basis=6-31G**;

        symmetry,nosym


        geometry=geom.xyz

        {{uhf;accu,5}}
        {{{method};{cas}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;}}

        show, energy
        table, energy
        save,enr.res,new
        {{table,____; noprint,heading}}
        '''.format(method=eInfo['method'],
                    state=eInfo['state'],
                    wf =  eInfo['wf'],
                    cas = eInfo['cas']))


    state = int(eInfo['state'])
    ########################
    #   Followinf is only for analytical NACT, 
    if nInfo['method'] == 'analytical': 
        nInfo['nmethod'] = 'cpmcscf'

    nactPairs = [[i,j] for i in range(1,state+1) for j in range(i+1,state+1)]


    nactTemp= "\n\nbasis={}\n".format(eInfo['basis'])

    for ind in range(0,len(nactPairs),5):
        nactTemp += "\n{{{method};{cas}; wf,{wf};state,{state}; start,2140.2;\n".format(
                            method=eInfo['method'],
                            state=eInfo['state'],
                            wf =  eInfo['wf'],
                            cas = eInfo['cas'])

        pairChunk = nactPairs[ind:ind+5]
        for i,j in pairChunk:
            nactTemp += "{nmethod},nacm,{f}.1,{s}.1,accu=1.d-10,record=51{f}{s}.1;\n".format(
                        nmethod = nInfo['nmethod'],
                        f=i,s=j)
        nactTemp += "}\n"



    for i,j in nactPairs:
        nactTemp +=textwrap.dedent("""
            force;nacm,51{f}{s}.1;varsav
            table,gradx,grady,gradz;
            save,ananac{f}{s}.res,new;
            {{table,____;  noprint,heading}}

            """.format(f=i,s=j))

    molproGridTemplate += nactTemp + '\n---\n'

    with open('grid.com','w') as f:
        f.write(molproGridTemplate)

    return state, nactPairs


def parseData(geomFile, freqFile, wilsonFile):
    geomData = np.loadtxt(geomFile, 
        dtype={'names': ('names', 'mass', 'x','y','z'),'formats': ('S1', 'f8', 'f8','f8','f8')})
    atomNames= geomData['names']
    atomsMass= geomData['mass']
    equiGeom = np.array([geomData[i] for i in ['x','y','z']]).T

    freq = np.loadtxt(freqFile)
    wilson = np.loadtxt(wilsonFile)
    return [atomNames, atomsMass, equiGeom, freq, wilson]


def createGeometry(nAtoms, atomNames, geomData, msg='', outFile = 'geom.xyz',):
    tmp = " {}\n".format(nAtoms)
    tmp+= "Geometry file created from ADT program. %s \n"%msg
    for i in range(nAtoms):
        tmp += "{},{},{},{}\n".format(atomNames[i], *geomData[i])

    with open(outFile,"w") as f:
        f.write(tmp)


def parseResult(file):
    with open(file,"r") as f:
        dat = f.read().replace("D","E").strip().split("\n")[1:]
    dat = [map(float,i.strip().split()) for i in dat]
    return np.array(dat)





def removeUnwantedFiles():
    filesToRemove = ['equi.out', 'equi.xml', 'grid.out', 'grid.xml']
    for file in filesToRemove:
        os.remove(file)


def rotationMatrix(phi):
    return np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])