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
    gInfo['firstgrid'] = map(float, gInfo['firstgrid'].split(','))
    gInfo['secondgrid'] = map(float, gInfo['secondgrid'].split(','))
    return [eInfo, nInfo,vModes, gInfo]



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








def parseResult(file):
    with open(file,"r") as f:
        dat = f.read().replace("D","E").strip().split("\n")[1:]
    dat = [map(float,i.strip().split()) for i in dat]
    return np.array(dat)





def removeUnwantedFiles():
    filesToRemove = ['equi.out', 'equi.xml', 'grid.out', 'grid.xml']
    for file in filesToRemove:
        os.remove(file)



def writeFile(data, file):
    pass
