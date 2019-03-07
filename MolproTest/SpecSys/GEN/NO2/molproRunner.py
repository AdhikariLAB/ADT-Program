from molproFuncs import *
import subprocess,os,shutil


# Molpro scratch directory
scrdir = '/tmp/koushik'


eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
state, nactPairs = createTemplate(eInfo, nInfo)
atomNames, atomsMass, equiGeom, freq, wilson = parseData('equiGeom.dat', 'freq.dat', 'wilson.dat')

equiGeom = equiGeom*0.529177


nAtoms = len(atomNames)
tmpWil = wilson
nModes = len(freq)
wilson = wilson.reshape(nAtoms, 3, nModes)

freqInv = np.sqrt(hbar/(freq*cInvToTauInv))
massInv = np.sqrt(1/atomsMass)


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


# energyResult = 
# nactResult   = 
mas = np.array([14.006700]*3 + [15.999400]*9)




for phi in phi_grid[:1]:
    phiRad = np.deg2rad(phi)

    for rho in rho_grid[:2]:
        shutil.copy('molpro.wfu',scrdir)

        qCord  = np.zeros(nModes)
        qCord[vModes[0]] = rho*np.cos(phiRad)
        qCord[vModes[1]] = rho*np.sin(phiRad)


        dsGeom  = np.einsum('ijk,k,k,i->ij',wilson, freqInv, qCord, massInv)
        curGeom = equiGeom+dsGeom

        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        createGeometry(nAtoms, atomNames, curGeom, msg)

        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==0:
            'Job successful '+msg
        else : 
            'Job unsuccessful'+msg 
            continue

        enrData = parseResult('enr.res').flatten()
        print rho, phi, enrData
        for i,j in nactPairs: 
            file = 'ananac{}{}.res'.format(i,j)
            nacts = parseResult(file).flatten()



            # nacts = np.abs(nacts)

            mtau = angtobohr*nacts*np.sqrt(hbar/mas)
            
            tauq1 = np.sum(mtau*tmpWil[:,4]/np.sqrt(freq[4]))
            tauq2 = np.sum(mtau*tmpWil[:,5]/np.sqrt(freq[5]))

            taurho = tauq1*np.cos(phiRad) + tauq2*np.sin(phiRad)

            # tau = np.einsum('ijk,ij,i,k->k',wilson,nacts,massInv,freqInv)

            print rho, phi,i,j, taurho



