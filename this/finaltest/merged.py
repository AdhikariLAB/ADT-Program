import os
import sys
import shutil
import textwrap
import subprocess
import numpy as np
import ConfigParser
from glob import glob
from datetime import datetime


#DON'T forget to remove all `noexec` from the ddr template


class Base():
    '''
        A base class containing the common methods to be used in both cases of spectroscopic and scattering
    '''
    angtobohr = 1.8897259886
    hbar = 0.063508
    cInvToTauInv = 0.001883651

    def sin(self, x):
        """ A sin function that directly takes degree as input unlike numpy"""
        return np.sin(np.deg2rad(x))

    def cos(self, x):
        """ A cos function that directly takes degree as input unlike numpy"""
        return np.cos(np.deg2rad(x))

    def createAnaTemplate(self):
        ''''Creates the molpro template files analytical job'''
        molproTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job.
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym


            geometry=geom.xyz

            {{{method};{cas}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;}}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}
            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        method=self.eInfo['method'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))


        nactTemp= "\n\nbasis={}\n".format(self.nInfo['basis'])

        for ind in range(0,len(self.nactPairs),5):
            nactTemp += "\n{{mcscf;{cas}; wf,{wf};state,{state}; start,2140.2;\n".format(
                                state=self.eInfo['state'],
                                wf =  self.eInfo['wf'],
                                cas = self.eInfo['cas'])

            pairChunk = self.nactPairs[ind:ind+5]
            forceTemp = ''
            for count,pair in enumerate(pairChunk, start=1):
                f,s = pair
                nactTemp += "cpmcscf,nacm,{f}.1,{s}.1,record=51{n:02}.1;\n".format(
                           f=f, s=s, n=count)

                forceTemp +=textwrap.dedent("""
                force;nacm,51{n:02}.1;varsav
                table,gradx,grady,gradz;
                save,ananac{f}{s}.res,new;
                {{table,____;  noprint,heading}}

                """.format(f=f, s=s, n=count))
            nactTemp += "}\n" + forceTemp


        molproTemplate += nactTemp + '\n---\n'

        with open('grid.com','w') as f:
            f.write(molproTemplate)



        molproInitTemplate = textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job
            memory,{memory}
            file,2,molpro_init.wfu,new;

            basis={basis};

            symmetry,nosym


            geometry=geom.xyz
            {{uhf}}
            {{{method};{cas}; wf,{wf};state,{state};orbital,2140.2;}}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}
            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        method=self.eInfo['method'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))

        with open('init.com', 'w') as f:
            f.write(molproInitTemplate)

    def createDdrTemplate(self):
        ''''Creates the molpro template files ddr job'''
        molproTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job for analytical job.
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom1.xyz

            {{multi;{cas}; wf,{wf};state,{state};orbital,2140.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;noexc}}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}

            !for +dr
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;noexc}}
            {{ci;trans,6000.2,6001.2;dm,8001.2}}


            !for -dr
            symmetry,nosym
            geometry=geom3.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;noexc}}
            {{ci;trans,6000.2,6002.2;dm,8002.2}}



            !for +dp
            symmetry,nosym
            geometry=geom4.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6003.2;noexc}}
            {{ci;trans,6000.2,6003.2;dm,8003.2}}


            !for -dp
            symmetry,nosym
            geometry=geom5.xyz
            {{multi;{cas}; wf,{wf};state,{state};start,2140.2;orbital,2244.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6004.2;noexc}}
            {{ci;trans,6000.2,6004.2;dm,8004.2}}


            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))


        nactTemp= ''

        for i,j in self.nactPairs:
            nactTemp+=textwrap.dedent(''' 
                !for taur     
                {{ddr,{dt},2140.2,2241.2,8001.2;state, {j}.1,{i}.1}}
                nacmepv=nacme

                {{ddr,-{dt},2140.2,2242.2,8002.2;state, {j}.1,{i}.1}}
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

                
                '''.format(dt=self.dt,dp=self.dp,i=i,j=j))



        molproTemplate += nactTemp + '\n---\n'

        with open('grid.com','w') as f:
            f.write(molproTemplate)

        molproInitTemplate = textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job.
            memory,{memory}
            file,2,molpro_init.wfu,new;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom.xyz

            {{uhf}}
            {{multi;{cas}; wf,{wf};state,{state};orbital,2140.2;}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;noexc}}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas = self.eInfo['cas']))



        with open('init.com', 'w') as f:
            f.write(molproInitTemplate)

    def makeGrid(self, ll):
        return np.arange(ll[0], ll[1]+ll[2], ll[2])

    def parseData(self, atomFile):
        ''' Parse atom names and masses from data file'''
        atomData = np.loadtxt(atomFile, 
            dtype={'names': ('names', 'mass'),'formats': ('S1', 'f8')})
        self.atomNames = atomData['names']
        self.atomMass  = atomData['mass']

    def parseConfig(self, conFigFile):
        ''' 
        Parses configuration for running molpro from the molpro config file provided
        and sets up different methods and attributes relavant to the configuration 
        '''
        scf = ConfigParser.SafeConfigParser()
        scf.read(conFigFile)

        molInfo =  dict(scf.items('molInfo'))
        self.scrdir = molInfo['scrdir']
        self.memory = molInfo['memory']

        self.eInfo = dict(scf.items('eInfo'))
        self.nInfo = dict(scf.items('nInfo'))
        gInfo = dict(scf.items('gInfo'))

        self.state = int(self.eInfo['state'])
        self.nactPairs = [[i,j] for i in range(1,self.state+1) for j in range(i+1,self.state+1)]
        self.nTau = len(self.nactPairs)
        #or should I just take as rho/theta/phi grid?
        self.rhoList = map(float, gInfo['rho'].split(','))
        self.thetaList = map(float, gInfo['theta'].split(','))
        self.phiList = map(float, gInfo['phi'].split(','))
        self.phiGrid = self.makeGrid(self.phiList)

        spec = self.__class__.__name__=='Spectroscopic'
        if spec:
            mInfo = dict(scf.items('mInfo'))
            self.vModes = [int(i)-1 for i in mInfo['varying'].split(',')]
            if len(self.rhoList) == 1:
                raise Exception('Give a grid rho value for spectroscopic calculation')
            self.rhoGrid = self.makeGrid(self.rhoList)


        else:
            if len(self.rhoList) != 1:
                raise Exception('Give a fixed rho value for scattering calculation')
            self.rho = self.rhoList[0]
            self.thetaGrid = self.makeGrid(self.thetaList)


        if self.nInfo['method']=='cpmcscf':
            self.createAnaTemplate()
            self.getTau = self.getTauAna
            self.createGridGeom = self.createOneGeom


        elif self.nInfo['method']=='ddr':
            if spec:
                self.dr = self.dt = float(gInfo['dr']) # in case of spectroscopic dt means dr
            else:
                self.dt = float(gInfo['dt'])
            self.dp = float(gInfo['dp'])
            self.createDdrTemplate()
            self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr

        else:
            sys.exit('%s Not a proper option'%self.nInfo['method'])
        self.logFile = open('adt_molpro.log', 'w', buffering=0)

    def parseResult(self, file):
        ''' Parses result from the ouptpu .res files'''
        with open(file,"r") as f:
            dat = f.read().replace("D","E").strip().split("\n")[1:]
        dat = [map(float,i.strip().split()) for i in dat]
        return np.array(dat)

    def writeFile(self, file, data):
        ''' Writes output data in plain txt'''
        file = open(file,'wb')
        for tp in np.unique(data[:,0]):
            np.savetxt( file, data[data[:,0]==tp] ,delimiter="\t", fmt="%.8f")
            file.write("\n")

    def interp(self, file ):
        ''' Fills the missing values in output file using a 1D interpolation '''
        data = np.loadtxt(file)
        # grid2 = self.phiGrid
        grid1 = np.unique(data[:,0])   # numpy precision error
        grid2 = np.unique(data[:,1])

        lim = data.shape[1]  # columns in data
        res = np.array([], dtype=np.float64).reshape(0,lim)

        for g1 in grid1:

            g1Block = data[data[:,0] == g1]

            res1 = np.column_stack([np.full(grid2.shape, g1), grid2])
            gx = g1Block[:,1]
            for col in range(2, lim):
                gy = g1Block[:,col]
                gny = np.interp(grid2, gx, gy)
                res1 = np.column_stack([res1, gny])
            res = np.vstack([res, res1])



        return res

    def removeFiles(self, allOut = False):
        ''' Removes unwanted files after each time running molpro'''
        # allOut True means include all the previous out files
        # by default it just include the current out files
        patterns = ['*.res', '*.xyz', '*.xml*', '*.out']
        files = []
        if allOut: patterns += ['*.wfu', '*.com']
        for i in patterns : 
            files.extend(glob(i))
        for file in files : 
            os.remove(file)

    def msg(self, m, cont=False):
        ''' Writes info in the log files'''
        if not cont : 
            m= "{:<30}{}".format(datetime.now().strftime("[%d-%m-%Y %I:%M:%S %p]"), m)
        else:
            m+='\n'
        self.logFile.write(m)

    def incompleteJobs(self, gridn1, g1, gridn2, g2):
        '''Saves the gemoetry and out file in a different folder '''
        # do this step previouslu
        path = '{}_{}_{}_{}'.format(gridn1, g1, gridn2, g2)
        path = os.path.join('IncompleteJobs', path)
        os.mkdir(path)
        files = glob('*.xyz')+ glob('*.out')
        for file in files:
            shutil.move(file, path)
        for file in glob('*.xml'): os.remove(file)


    def runThisMolpro(self, grid1, gridn1, grid2, gridn2, filEe, fileN1, fileN2):
        '''Runs the molpro for each gridpoints '''
        # subprocess calls blocks system I/O buffer, so the I/Os (sometimes) have to be flushed out explicitely 
        # open files to store result
        filee  = open(filEe, 'w', buffering=1)
        filen1 = open(fileN1,'w', buffering=1)
        filen2 = open(fileN2,'w', buffering=1)

        self.equiRun()
        self.removeFiles()
        #grid1 is theta or rho
        #grid2 is phi
        #create folder for incomplete jobs , delete if already exists
        if os.path.isdir('IncompleteJobs'): shutil.rmtree('IncompleteJobs')
        os.mkdir('IncompleteJobs')

        done = False # why am I using this?

        for g2 in grid2:
            shutil.copy('molpro_init.wfu', 'molpro.wfu')      # copy the wfu from the initial/equilibrium job
            for g1 in grid1: 

                #now if rho 0 is provided run the rho=0=phi step only once
                if (g1==0) & (g2!=0) & (not done):
                    continue # don't run the job for other rhos
                else:
                    done = True
                
                self.createGridGeom(g1, g2)
                self.msg( 'Running molpro job for {} = {}, {} = {}.......'.format(gridn1, g1, gridn2, g2))
                shutil.copy('molpro.wfu',self.scrdir)
                exitcode = subprocess.call(['molpro',"-d", self.scrdir,'-W .','grid.com'])
                if exitcode==0:
                    self.msg( 'Job successful.', cont=True)
                else:
                    self.msg( 'Job failed.', cont=True)
                    self.incompleteJobs(gridn1, g1, gridn2, g2)
                    continue

                enrData = self.parseResult('enr.res').flatten()
                tau1, tau2 = self.getTau(g1, g2)
                np.savetxt(filee,  np.append([g1,g2],enrData)[None], fmt='%.8f', delimiter='\t')
                np.savetxt(filen1, np.append([g1,g2],tau1)[None],  fmt='%.8f', delimiter='\t')
                np.savetxt(filen2, np.append([g1,g2],tau2)[None],  fmt='%.8f', delimiter='\t')
                self.removeFiles()
            filee.write('\n')
            filen1.write('\n')
            filen2.write('\n')
        self.removeFiles(allOut=True)
        self.msg('All molpro jobs done.')

        scat = self.__class__.__name__=='Scattering'
        #fill the missing values by 1D interpolation
        # '_mod' suffix means files with filled data
        for file in [filEe, fileN1, fileN2]:
            dat = self.interp(file)
            dat[:,1] = np.deg2rad(dat[:,1])
            if scat: dat[:,0] = np.deg2rad(dat[:,0])  # for scattering also transform the column 0
            if file == filEe: dat[:,2:] -= np.loadtxt('equienr.dat')[0]
            self.writeFile(file.replace('.dat', '_mod.dat'), dat)



class Spectroscopic(Base):
    ''' Inherited from the Base class, this class containg necessary methods
     for running molpro for a Spectroscopic system'''
    def __init__(self, conFigFile ,atomFile, geomFile , freqFile , wilsonFile ):
        self.parseData(atomFile)
        self.parseSData(geomFile, freqFile, wilsonFile)
        self.parseConfig(conFigFile)

    def parseSData(self, geomFile, freqFile, wilsonFile):
        '''Parses equilibrium geometry, frequencies and the wilson matrix data for a sceptroscopic system'''
        self.equiGeom = np.loadtxt(geomFile)
        freq = np.loadtxt(freqFile)
        wilson = np.loadtxt(wilsonFile)

        wilson = wilson.reshape(self.equiGeom.shape[0], 3, freq.shape[0])
        freqInv = np.sqrt(self.hbar/(freq*self.cInvToTauInv))
        massInv = np.sqrt(1/self.atomMass)
        self.wilFM = np.einsum('ijk,k,i->ijk',wilson,freqInv,massInv)

    def createOneGeom(self, rho, phi, outFile = 'geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given rho , phi'''
        nModes = self.wilFM.shape[2]
        qCord  = np.zeros(nModes)
        qCord[self.vModes[0]] = rho*self.cos(phi)
        qCord[self.vModes[1]] = rho*self.sin(phi)

        curGeom  = self.equiGeom+np.einsum('ijk,k->ij',self.wilFM, qCord)
        msg = 'for Rho = {}, Phi = {}'.format(rho, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])

        with open(outFile,"w") as f:
            f.write(tmp)

    def createAllGeom(self, rho, phi):
        ''' Creates 5 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom(rho,  phi,  'geom1.xyz')
        self.createOneGeom(rho+self.dr, phi, 'geom2.xyz')
        self.createOneGeom(rho-self.dr, phi, 'geom3.xyz')
        self.createOneGeom(rho, phi+self.dp,'geom4.xyz')
        self.createOneGeom(rho, phi-self.dp,'geom5.xyz')
    
    def parseNact(self, i,j,rho,phi):
        '''Calculates NACT Rho Phi from the output gradients '''
        file = 'ananac{}{}.res'.format(i,j)

        grads = self.angtobohr*self.parseResult(file)
        rotoMat = np.array([[self.cos(phi), self.sin(phi)], [-rho*self.sin(phi), rho*self.cos(phi)]])
        tau = np.einsum('ijk,ij,lk->l',self.wilFM[...,self.vModes], grads, rotoMat)
        return np.abs(tau)

    def getTauAna(self, rho, phi):
        '''Used in Analytical NACT calculation'''
        return np.stack([self.parseNact(i, j, rho, phi) 
                                    for i,j in self.nactPairs]).T

    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        return np.stack([self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs]).T

    def equiRun(self):
        '''Runs molpro for the equilibrium geometry'''
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program for equilibrium geometry.\n"
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *self.equiGeom[i])

        with open('geom.xyz',"w") as f:
            f.write(tmp)

        self.msg( "Running molpro job for equilibrium point.......")
        exitcode = subprocess.call(['molpro',"-d", self.scrdir ,'-W .','init.com'])
        if exitcode==0: 
            self.msg( 'Job successful', cont=True)
        else:
            self.msg( 'Job failed', cont=True)
            sys.exit('Molpro failed in equilibrium step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt='%.8f')
        self.removeFiles()



    def runMolpro(self):
        self.runThisMolpro(
            self.rhoGrid, 
            'Rho', 
            self.phiGrid, 
            'Phi', 
            'energy.dat', 
            'tau_rho.dat', 
            'tau_phi.dat' )



class Scattering(Base):
    ''' Inherited from the Base class, this class containg necessary methods
     for running molpro for a Scattering system'''
    def __init__(self, conFigFile, atomFile ):
        self.parseData(atomFile)
        self.parseConfig(conFigFile)

    def AreaTriangle(self,a,b,c):
        """ area of a tringle with sides a,b,c """
        ps = (a+b+c)/2.0
        ar = ps*(ps-a)*(ps-b)*(ps-c)
        # negative area due to round off errors set to zero
        if ar < 0.0:
            ar = 0.0
        ar = np.sqrt(ar)
        return ar

    def toJacobi(self,theta,phi):

        """ returns jacobi coordinates """
        m1, m2, m3 = self.atomMass

        M = m1 + m2 + m3
        mu = np.sqrt(m1*m2*m3/M)
        d1 = np.sqrt(m1*(m2+m3)/(mu*M))
        d2 = np.sqrt(m2*(m3+m1)/(mu*M))
        d3 = np.sqrt(m3*(m1+m2)/(mu*M))
        eps3 = 2 * np.arctan(m2/mu)
        eps2 = 2 * np.arctan(m3/mu)

        R1 = (1.0/np.sqrt(2.0))*self.rho*d3*np.sqrt(1.0+ self.sin(theta)*self.cos(phi+eps3)) # F-H2 distance
        R2 = (1.0/np.sqrt(2.0))*self.rho*d1*np.sqrt(1.0+ self.sin(theta)*self.cos(phi))      # H1-H2 distance
        R3 = (1.0/np.sqrt(2.0))*self.rho*d2*np.sqrt(1.0+ self.sin(theta)*self.cos(phi-eps2)) # F-H1 distance

        if R1 < 1e-10:
            R1 = 0.0
        if R2 < 1e-10:
            R2 = 0.0
        if R3 < 1e-10:
            R3 = 0.0

        area = self.AreaTriangle(R1,R2,R3)
        x = R2*R2 + R3*R3 - R1*R1
        y = 4.0*area
        Ang123 = np.arctan2(y,x)
        x2 = (0.0,0.0)
        x3 = (R2,0.0)
        x1 = (R3*np.cos(Ang123),R3*np.sin(Ang123))
        # these are non-mass scaled jacobi coords
        # r : (x3-x2)
        # R : x1 - com(x3,x2)
        # gamma : angle between r and R
        r = (R2,0.0)
        R = (R3*np.cos(Ang123) - m3*R2/(m2+m3) , R3*np.sin(Ang123))
        rs = np.linalg.norm(r)
        rc = np.linalg.norm(R)
        if rc < 1e-10:
            rc = 0.0

        rtil = (R2*m2/(m2+m3),0.0)
        drtil = np.linalg.norm(rtil)
        Areasmall = self.AreaTriangle(drtil,R1,rc)
        y = 4.0 * Areasmall
        x = drtil*drtil + rc*rc - R1*R1
        if np.fabs(x) < 1e-10:
            x = 0.0
        gamma = np.arctan2(y,x)
        return (rs, rc, gamma)

    def hyperToCart(self, theta, phi):
        rs, rc, gamma = self.toJacobi(theta, phi)
        p1 = [0, rc*np.cos(gamma),rc*np.sin(gamma)]
        p2 = [0, -rs/2.0 , 0.0 ]
        p3 = [0, rs/2.0  , 0.0 ]
        return np.array([p1,p2,p3])

    def createOneGeom(self, theta, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given theta , phi'''
        curGeom  = self.hyperToCart(theta, phi)
        msg = 'for Rho = {}, Theta={}, Phi = {}'.format(self.rho, theta, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])
        with open(outFile,"w") as f:
            f.write(tmp)

    def createAllGeom(self, theta, phi):
        ''' Creates 5 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom(theta,    phi,    'geom1.xyz')
        self.createOneGeom(theta+self.dt, phi,   'geom2.xyz')
        self.createOneGeom(theta-self.dt, phi,   'geom3.xyz')
        self.createOneGeom(theta,    phi+self.dp,'geom4.xyz')
        self.createOneGeom(theta,    phi-self.dp,'geom5.xyz')

    def parseNact(self, i,j,gradTheta, gradPhi):
        '''Calculates NACT Rho Phi from the output gradients '''
        file = 'ananac{}{}.res'.format(i,j)
        grads = self.angtobohr*self.parseResult(file)

        tauTheta = np.einsum('ij,ij', grads, gradTheta)
        tauPhi   = np.einsum('ij,ij', grads, gradPhi)
        return np.abs(np.array([tauTheta, tauPhi]))

    def getTauAna(self, theta, phi):
        '''Used in Analytical NACT calculation'''
        # What will be this values
        dTheta = self.thetaList[2]/100
        dPhi = self.phiList[2]/100

        gradThetaPlus  = self.hyperToCart(theta+dTheta, phi)
        gradThetaMinus = self.hyperToCart(theta-dTheta, phi)
        gradTheta      = (gradThetaPlus - gradThetaMinus)/2*dTheta

        gradPhiPlus  = self.hyperToCart(theta+dTheta, phi)
        gradPhiMinus = self.hyperToCart(theta-dTheta, phi)
        gradPhi      = (gradPhiPlus - gradPhiMinus)/2*dPhi 
        return np.vstack([self.parseNact(i,j,gradTheta, gradPhi) 
                                    for i,j in self.nactPairs]).T

    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        return np.vstack([self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs]).T


    def equiRun(self):
        '''Runs molpro for a initial geometry, i.e. theta=phi=0'''
        self.createOneGeom(0, 0)
        self.msg( "Running molpro job for initial point....." )
        sys.stdout.flush()
        exitcode = subprocess.call(['molpro',"-d", self.scrdir,'-W .','init.com'])
        if exitcode==0: 
            self.msg( 'Job successful', cont= True)
        else:
            self.msg( 'Job failed', cont= True)
            sys.exit('Molpro failed in initital step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt='%.8f')
        self.removeFiles()



    def runMolpro(self):
        self.runThisMolpro(
            self.thetaGrid, 
            'Theta', 
            self.phiGrid, 
            'Phi', 
            'energy.dat', 
            'tau_theta.dat', 
            'tau_phi.dat' )



if __name__ == "__main__":
    # s = Scattering('./scatterAna/molpro.config','./scatterAna/atomfile.dat')


    s = Scattering('./ScatterDdr/molpro.config','./ScatterDdr/atomfile.dat')


    # s = Spectroscopic('./specAan/molpro.config','./specAan/atomfile.dat', './specAan/geomfile.dat', './specAan/frequency.dat', './specAan/wilson.dat')


    # s = Spectroscopic('./specDdr/molpro.config','./specDdr/atomfile.dat', './specDdr/geomfile.dat', './specDdr/frequency.dat', './specDdr/wilson.dat')


    s.runMolpro()