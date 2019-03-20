from molproFuncs import *
import subprocess,os,shutil


# Molpro scratch directory
scrdir = '/tmp/koushik'


eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
state, nactPairs   = createTemplate(eInfo, nInfo)
atomNames, equiGeom, wilFM = parseData('equiGeom.dat', 'freq.dat', 'wilson.dat')

equiGeom = equiGeom*0.529177

rho_grid = np.arange(0.1,5.1, 0.1)
phi_grid = np.arange(1,181, 2)


nModes = wilFM.shape[2]
nTau = len(nactPairs)

energyResult  = np.array([], dtype=np.float64).reshape(0,state+2) 
nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
#assert atom number and mode relation




#################### Equilibrium step sta ################
createGeometry(atomNames, equiGeom, 'For equilibrium')


exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','equi.com'])
if exitcode==1: sys.exit('Molpro failed in equilibrium step')

equiData = parseResult('equienr.res')
print equiData.flatten()
#####################################################





def parseNact(i,j,rho,phi):
    file = 'ananac{}{}.res'.format(i,j)

    nacts = angtobohr*parseResult(file)
    rotoMat = np.array([[np.cos(phi), np.sin(phi)], [-rho*np.sin(phi), rho*np.cos(phi)]])
    tau = np.einsum('ijk,ij,lk->l',wilFM[...,vModes],nacts,rotoMat)

    return np.abs(tau)




for phi in phi_grid[:1]:
    phiRad = np.deg2rad(phi)

    for rho in rho_grid[:2]:
        shutil.copy('molpro.wfu',scrdir)

        qCord  = np.zeros(nModes)
        qCord[vModes[0]] = rho*np.cos(phiRad)
        qCord[vModes[1]] = rho*np.sin(phiRad)

        dsGeom  = np.einsum('ijk,k->ij',wilFM,qCord)
        curGeom = equiGeom+dsGeom

        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        createGeometry(atomNames, curGeom, msg)

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