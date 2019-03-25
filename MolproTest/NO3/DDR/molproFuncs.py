import numpy as np
import ConfigParser ,textwrap, sys, os


angtobohr = 1.8897259886
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



def createTemplate(eInfo, nInfo,dr,dp):
    molproEquiTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        memory,500,m
        file,2,molpro.wfu,new;

        basis=6-31G**;

        symmetry,nosym


        geometry=geom.xyz

        {{uhf;orbital,2200.2}}
        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};orbital,2140.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6000.2;noexc}}
        show,  energy
        table, energy
        save,equienr.res,new
        {{table,____; noprint,heading}}

        ---

        '''.format( state=eInfo['state'],
                    wf =  eInfo['wf'],
                    occ= eInfo['occ'],
                    closed = eInfo['closed'],
                    core = eInfo['core']))

    with open('equi.com' ,'w') as f:
        f.write(molproEquiTemplate)



    molproGridTemplate=textwrap.dedent('''
        ***, Molpro template created from ADT program for "WhatEver"
        memory,500,m
        file,2,molpro.wfu;

        basis=6-31G**;

        symmetry,nosym

        geomtyp=xyz
        geometry=geom1.xyz

        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};orbital,2140.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6000.2;noexc}}

        show, energy
        table, energy
        save,enr.res,new
        {{table,____; noprint,heading}}

        !for +dr
        symmetry,nosym
        geometry=geom2.xyz
        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};start,2140.2;orbital,2241.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6001.2;noexc}}
        {{ci;trans,6000.2,6001.2;dm,8001.2}}


        !for -dr
        symmetry,nosym
        geometry=geom3.xyz
        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6002.2;noexc}}
        {{ci;trans,6000.2,6002.2;dm,8002.2}}



        !for +dp
        symmetry,nosym
        geometry=geom4.xyz
        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6003.2;noexc}}
        {{ci;trans,6000.2,6003.2;dm,8003.2}}


        !for -dp
        symmetry,nosym
        geometry=geom5.xyz
        {{multi;occ,{occ};closed,{closed}; wf,{wf};state,{state};start,2140.2;orbital,2244.2;}}
        {{mrci; occ,{occ};closed,{closed};core,{core}; wf,{wf};state,{state};save,6004.2;noexc}}
        {{ci;trans,6000.2,6004.2;dm,8004.2}}


        '''.format( state=eInfo['state'],
                    wf =  eInfo['wf'],
                    occ= eInfo['occ'],
                    closed = eInfo['closed'],
                    core = eInfo['core']))


    state = int(eInfo['state'])
    ########################
    #   Following is only for ddr NACT, 






    nactPairs = [[i,j] for i in range(1,state+1) for j in range(i+1,state+1)]

    nactTemp= ''

    for i,j in nactPairs:
        nactTemp+=textwrap.dedent(''' 
            !for taur     
            {{ddr,{dr},2140.2,2241.2,8001.2;state, {j}.1,{i}.1}}
            nacmepv=nacme

            {{ddr,-{dr},2140.2,2242.2,8002.2;state, {j}.1,{i}.1}}
            nacmemv=nacme

            nacmr = 0.5*(nacmepv+ nacmemv)

            !for taup
            {{ddr,{dp},2140.2,2243.2,8003.2;state, {j}.1,{i}.1}}
            nacmepv=nacme

            {{ddr,-{dp},2140.2,2244.2,8004.2;state, {j}.1,{i}.1}}
            nacmemv=nacme
            nacmp = 0.5*(nacmepv+ nacmemv)


            table, nacmr,nacmp
            save,ddrnact{i}{j}.res,new;

            
            '''.format(dr=dr,dp=dp,i=i,j=j))



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
    equiGeom = equiGeom

    freq = np.loadtxt(freqFile)
    wilson = np.loadtxt(wilsonFile)

    wilson = wilson.reshape(atomNames.shape[0], 3, freq.shape[0])
    freqInv = np.sqrt(hbar/(freq*cInvToTauInv))
    massInv = np.sqrt(1/atomsMass)
    wilFM = np.einsum('ijk,k,i->ijk',wilson,freqInv,massInv)
    return [atomNames, equiGeom, wilFM]





def createEquiGeometry(atomNames, geomData, msg=''):
    nAtoms = len(atomNames)
    tmp = " {}\n".format(nAtoms)
    tmp+= "Geometry file created from ADT program. %s \n"%msg
    for i in range(nAtoms):
        tmp += "{},{},{},{}\n".format(atomNames[i], *geomData[i])

    with open('geom.xyz',"w") as f:
        f.write(tmp)




def createAllGridGeom(atomNames, equiGeom,wilFM,vModes, rho, phi, dr, dp):
    createGridGeom(atomNames, equiGeom, wilFM,vModes, rho, phi,'geom1.xyz')
    createGridGeom(atomNames, equiGeom, wilFM,vModes, rho+dr, phi,'geom2.xyz')
    createGridGeom(atomNames, equiGeom, wilFM,vModes, rho-dr, phi,'geom3.xyz')
    createGridGeom(atomNames, equiGeom, wilFM,vModes, rho, phi+dp,'geom4.xyz')
    createGridGeom(atomNames, equiGeom, wilFM,vModes, rho, phi-dp,'geom5.xyz')


def createGridGeom(atomNames, equiGeom,wilFM,vModes, rho, phi,outFile):
    nModes = wilFM.shape[2]
    qCord  = np.zeros(nModes)
    qCord[vModes[0]] = rho*np.cos(phi)
    qCord[vModes[1]] = rho*np.sin(phi)

    curGeom  = equiGeom+np.einsum('ijk,k->ij',wilFM, qCord)

    nAtoms = len(atomNames)
    tmp = " {}\nGeometry file created from ADT program. Rho={} Phi={}\n".format(nAtoms,rho,phi)

    for i in range(nAtoms):
        tmp += "{},{},{},{}\n".format(atomNames[i], *curGeom[i])

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


