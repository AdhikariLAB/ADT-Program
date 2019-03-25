from molproFuncs import *
import subprocess,os,shutil

scrdir = '/tmp/koushik'  #take it from config

eInfo, nInfo, vModes, gInfo = parseConfig('molpro.config')
atomNames, equiGeom, wilFM  = parseData('equiGeom.dat', 'freq.dat', 'wilson.dat')



#check for type of nact
if nInfo['method']=='ddr':
    from ddrFuncs import *
    dr = gInfo['firstgrid'][2]/100              #setting the dr dp as 1/100 times the stepsize
    dp = gInfo['secondgrid'][2]/100
    parseNact = parseResult

elif nInfo['method']=='analytical':
    from anaFuncs import *
    createEquiGeom = createGeometry
    dr = None
    dp = None
    # define dr and dp as none
else :
    sys.exit()

state, nactPairs   = createTemplate(eInfo, nInfo,dr,dp)


# equiGeom = equiGeom*0.529177
#create this from config file
rho_grid = np.arange(*gInfo['firstgrid'])
phi_grid = np.arange(*gInfo['secondgrid'])
nModes = wilFM.shape[2]
nTau = len(nactPairs)



#assert atom number and mode relation
energyResult  = np.array([], dtype=np.float64).reshape(0,state+2) 
nactRhoResult = np.array([], dtype=np.float64).reshape(0,nTau+2)
nactPhiResult = np.array([], dtype=np.float64).reshape(0,nTau+2)





################### Equilibrium step ################
createEquiGeom(atomNames, equiGeom, 'For equilibrium')


exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','equi.com'])
if exitcode==1: sys.exit('Molpro failed in equilibrium step')

equiData = parseResult('equienr.res')
print equiData.flatten()
######################################################



### Warning !!!! Check for degree and radian confusion

for phi in phi_grid:

    for rho in rho_grid[:1]:

        createGridGeom(atomNames, equiGeom,wilFM, vModes, rho, phi, dr, dp)
        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)


        shutil.copy('molpro.wfu',scrdir)
        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==0:
            print 'Job successful  '+msg
        else : 
            print 'Job unsuccessful'+msg 
            continue


        enrData = parseResult('enr.res').flatten()
        #fix parsenact for analytical and ddr
        tauRho, tauPhi = np.stack((parseNact('ddrnact{}{}.res'.format(i,j)) for i,j in nactPairs)).T


        energyResult  = np.vstack((energyResult,  np.append([rho,phi],enrData)))
        nactRhoResult = np.vstack((nactRhoResult, np.append([rho,phi],tauRho)))
        nactPhiResult = np.vstack((nactPhiResult, np.append([rho,phi],tauPhi)))

print energyResult, nactRhoResult, nactPhiResult