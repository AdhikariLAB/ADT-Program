from molproFuncs import *
import subprocess,os,shutil


# Molpro scratch directory
scrdir = '/tmp/koushik'


eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
state, nactPairs   = createTemplate(eInfo, nInfo)
atomNames, atomsMass, equiGeom, freq, wilson = parseData('equiGeom.dat', 'freq.dat', 'wilson.dat')

equiGeom = equiGeom*0.529177


nAtoms = len(atomNames)
# tmpWil = wilson
nModes = len(freq)
nTau = len(nactPairs)
wilson = wilson.reshape(nAtoms, 3, nModes)
freqInv = np.sqrt(hbar/(freq*cInvToTauInv))
massInv = np.sqrt(1/atomsMass)
wilFreqMass = np.einsum('ijk,k,i->ijk',wilson,freqInv,massInv)

#assert atom number and mode relation




#################### Equilibrium step sta ################
createGeometry(nAtoms, atomNames, equiGeom, 'For equilibrium')


exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','equi.com'])
if exitcode==1: sys.exit('Molpro failed in equilibrium step')

equiData = parseResult('equienr.res')
print equiData.flatten()
#####################################################

rho_grid = np.arange(0.1,5.1, 0.1)
phi_grid = np.arange(1,181, 2)


energyResult  = np.array([], dtype=np.float64).reshape(0,state+2) 
nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)


def parseNact(i,j,rho,phi):
    file = 'ananac{}{}.res'.format(i,j)

    nacts = angtobohr*parseResult(file)
    tau = np.einsum('ijk,ij->k',wilFreqMass,nacts)

    rotoMat = np.array([[np.cos(phi), np.sin(phi)], [-rho*np.sin(phi), rho*np.cos(phi)]])
    tau = np.dot(rotoMat,tau[vModes])
    return np.abs(tau)




for phi in phi_grid[:1]:
    phiRad = np.deg2rad(phi)

    for rho in rho_grid[:2]:
        shutil.copy('molpro.wfu',scrdir)

        qCord  = np.zeros(nModes)
        qCord[vModes[0]] = rho*np.cos(phiRad)
        qCord[vModes[1]] = rho*np.sin(phiRad)


        dsGeom  = np.einsum('ijk,k->ij',wilFreqMass,qCord)
        curGeom = equiGeom+dsGeom

        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        createGeometry(nAtoms, atomNames, curGeom, msg)

        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==0:
            print 'Job successful  '+msg
        else : 
            print 'Job unsuccessful'+msg 
            continue

        enrData = parseResult('enr.res').flatten()
        tauRho, tauPhi = np.stack((parseNact(i,j,rho,phiRad) for i,j in nactPairs)).T

        energyResult  = np.vstack((energyResult,  np.append([rho,phiRad],enrData)))
        nactRhoResult = np.vstack((nactRhoResult, np.append([rho,phiRad],tauRho)))
        nactPhiResult = np.vstack((nactPhiResult, np.append([rho,phiRad],tauPhi)))

print energyResult, nactRhoResult, nactPhiResult