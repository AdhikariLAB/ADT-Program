from molproFuncs import *
import subprocess,os,shutil


# Molpro scratch directory
scrdir = '/tmp/koushik'

dp = 2/100.0
dr = 0.1/100.0
eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
state, nactPairs   = createTemplate(eInfo, nInfo,dr,dp)
atomNames, equiGeom, wilFM = parseData('equiGeom.dat', 'freq.dat', 'wilson.dat')

equiGeom = equiGeom*0.529177

rho_grid = np.arange(0.1,5.1, 0.1)
phi_grid = np.arange(1,181, 2)

# tmpWil = wilson
nModes = wilFM.shape[2]
nTau = len(nactPairs)

energyResult  = np.array([], dtype=np.float64).reshape(0,state+2) 
nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
# #assert atom number and mode relation




# #################### Equilibrium step ################
createEquiGeometry(atomNames, equiGeom, 'For equilibrium')


exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','equi.com'])
if exitcode==1: sys.exit('Molpro failed in equilibrium step')

equiData = parseResult('equienr.res')
print equiData.flatten()
# #####################################################









# dp must be in radian
dp = np.deg2rad(dp)

for phi in phi_grid[:1]:
    phiRad = np.deg2rad(phi)

    for rho in rho_grid[:1]:
        shutil.copy('molpro.wfu',scrdir)

        createAllGridGeom(atomNames, equiGeom,wilFM, vModes, rho, phiRad, dr, dp)
        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==0:
            print 'Job successful  '+msg
        else : 
            print 'Job unsuccessful'+msg 
            continue

        enrData = parseResult('enr.res').flatten()
        tauRho, tauPhi = np.stack((parseResult('ddrnact{}{}.res'.format(i,j)) for i,j in nactPairs)).T
        print tauRho
        print tauPhi

        energyResult  = np.vstack((energyResult,  np.append([rho,phiRad],enrData)))
        nactRhoResult = np.vstack((nactRhoResult, np.append([rho,phiRad],tauRho)))
        nactPhiResult = np.vstack((nactPhiResult, np.append([rho,phiRad],tauPhi)))

# print energyResult, nactRhoResult, nactPhiResult