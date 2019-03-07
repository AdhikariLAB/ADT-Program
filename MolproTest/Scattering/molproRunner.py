from molproFuncs import *
import subprocess,os,shutil


# Molpro scratch directory
scrdir = '/tmp/koushik'


eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
state, nactPairs   = createTemplate(eInfo, nInfo)


# equiGeom = equiGeom*0.529177




nTau = len(nactPairs)

#assert atom number and mode relation



rho = 
pStart, pEnd, pStep = 
tStart, tEnd, tStep = 



#####################################################

theta_grid = np.arange(pStart, pEnd+ pStep, pStep )
phi_grid = np.arange(tStart, tEnd + tStep, tStep )


energyResult  = np.array([], dtype=np.float64).reshape(0,state+2) 
nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)


for theta in theta_grid : 

    for phi in phi_grid:
        shutil.copy('molpro.wfu',scrdir)   #1st step handling 
        curGeom = hyperToCart(rho, theta, phi)
        msg = 'For Rho = {}, Theta = {}, Phi = {}'.format(rho, theta, phi)
        createGeometry(atomNames, curGeom, msg)
        
        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==0:
            print 'Job successful  '+msg
        else : 
            print 'Job unsuccessful'+msg 
            continue




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