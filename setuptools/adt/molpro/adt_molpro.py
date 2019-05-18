
__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''
import re
import os
import sys
import shutil
import textwrap
import subprocess
import numpy as np
import ConfigParser
from glob import glob
from datetime import datetime
from adtmod import adt



class Base():
    '''
        A base class containing the common methods to be used in both cases of spectroscopic and scattering
    '''
    angtobohr    = 1.8897259886
    hbar         = 0.063508
    cInvToTauInv = 0.001883651
    bohrtoang    = 0.529177

    def sin(self, x):
        """ A sin function that directly takes degree as input unlike numpy"""
        return np.sin(np.deg2rad(x))

    def cos(self, x):
        """ A cos function that directly takes degree as input unlike numpy"""
        return np.cos(np.deg2rad(x))

    def interpolate(self,x,y,newx):
        diff = adt.spline(x,y)
        return np.array([ adt.splint(x,y,diff,xi) for xi in newx])

    def createAnaTemplate(self):
        ''''Creates the molpro template files analytical job'''
        cas = self.eInfo['cas']

        #remove core from the cas
        try:
            core = re.search( '[0-9a-zA-Z,;](core,\d+;*)', cas).group(1)
            cas = cas.replace(core, '').strip(';')
        except:
            pass


        energyLine ='{{mcscf;{cas}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra}}}'.format(
                        state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas=cas,
                        extra= self.eInfo['multi_extra'])

        if self.eInfo['method'] == 'mrci':
            energyLine+='''
            {{mrci;{cas}; wf,{wf};state,{state};{extra}}}
            '''.format(state=self.eInfo['state'],
                        wf =  self.eInfo['wf'],
                        cas=self.eInfo['cas'],
                        extra = self.eInfo['mrci_extra'])

        molproTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job.
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym


            geometry=geom.xyz
            {enrl}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}
            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        enrl = energyLine))


        nactTemp= "\n\nbasis={}\n".format(self.nInfo['basis'])

        for ind in range(0,len(self.nactPairs),5):
            nactTemp += "\n{{mcscf;{cas}; wf,{wf};state,{state}; start,2140.2;{extra};\n".format(
                                state=self.eInfo['state'],
                                wf =  self.eInfo['wf'],
                                cas=cas,
                                extra= self.eInfo['multi_extra'])

            pairChunk = self.nactPairs[ind:ind+5]
            forceTemp = ''
            for count,pair in enumerate(pairChunk, start=1):
                f,s = pair
                nactTemp += "cpmcscf,nacm,{f}.1,{s}.1,record=51{n:02}.1;{extra};\n".format(
                    f=f, s=s, n=count, extra=self.nInfo['nact_extra'])

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
            {{hf}}
            {enrl}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            ---

            '''.format(memory = self.memory,
                        basis = self.eInfo['basis'],
                        enrl=energyLine,
                        hf = self.eInfo['hf']))

        with open('init.com', 'w') as f:
            f.write(molproInitTemplate)

    def createDdrTemplate(self):
        ''''Creates the molpro template files ddr job'''

        cas = self.eInfo['cas']

        #remove core from the cas
        try:
            core = re.search( '[0-9a-zA-Z,;](core,\d+;*)', cas).group(1)
            cas = cas.replace(core, '').strip(';')
        except:
            pass

        #mrci has to be done for ddr nact calculation
        energyLine = """
            {{mcscf;{cas1}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;{extra2}}}
            """.format(state   =self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'])


        molproTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job for analytical job.
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom1.xyz

            {enrl}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}

            !for +dr
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas1} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;{extra2}}}
            {{ci;trans,6000.2,6001.2;dm,8001.2;{extra3}}}


            !for -dr
            symmetry,nosym
            geometry=geom3.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;{extra2}}}
            {{ci;trans,6000.2,6002.2;dm,8002.2;{extra3}}}



            !for +dp
            symmetry,nosym
            geometry=geom4.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6003.2;{extra2}}}
            {{ci;trans,6000.2,6003.2;dm,8003.2;{extra3}}}


            !for -dp
            symmetry,nosym
            geometry=geom5.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2244.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6004.2;{extra2}}}
            {{ci;trans,6000.2,6004.2;dm,8004.2;{extra3}}}


            '''.format( memory = self.memory,
                        basis  = self.eInfo['basis'],
                        enrl   = energyLine,
                        state  = self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'],
                        extra3 = self.nInfo['nact_extra']))


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
            ***, Molpro template created from ADT program for ddr job.
            memory,{memory}
            file,2,molpro_init.wfu,new;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom.xyz

            {{hf}}
            {enrl}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        enrl=energyLine,
                        hf = self.eInfo['hf']))



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

    def parseConfig(self, scf):
        ''' 
        Parses configuration for running molpro from the provided molpro config file 
        and sets up different methods and attributes relavant to the configuration 
        '''


        spec = self.__class__.__name__ == 'Spectroscopic'

        molInfo =  dict(scf.items('molInfo'))
        self.scrdir = molInfo['scrdir']
        self.memory = molInfo['memory']
        try:
            self.proc = molInfo['processor']
        except KeyError:
            self.proc = '1'

        self.eInfo = dict(scf.items('eInfo'))
        self.nInfo = dict(scf.items('nInfo'))

        #set some default values for the optional arguments

        if not 'multi_extra' in self.eInfo.keys():
            self.eInfo['multi_extra'] = ''
        if not 'mrci_extra' in self.eInfo.keys():
            self.eInfo['mrci_extra'] = ''
        if not 'nact_extra' in self.nInfo.keys():
            self.nInfo['nact_extra'] = ''


        # set up initial HF calculation according to 'restricted' parameter
        self.eInfo['hf'] = 'uhf'
        if self.eInfo['restricted'].lower() == 'true':
            self.eInfo['hf'] = 'hf'
        if 'uhf_extra' in self.eInfo.keys():
            self.eInfo['hf'] += ';' + self.eInfo['uhf_extra']



        gInfo = dict(scf.items('gInfo'))

        self.state = int(self.eInfo['state'])
        self.nactPairs = [[i, j] for j in range(2, self.state+1) for i in range(1, j)]

        # self.nactPairs = [[i,j] for i in range(1,self.state+1) for j in range(i+1,self.state+1)]
        self.nTau = len(self.nactPairs)

        self.phiList = map(float, gInfo['phi'].split(','))
        self.phiGrid = self.makeGrid(self.phiList)

        # there will be there types 
        # i. normal mode/ spectroscopic
        # ii. scattering hyperspherical
        # iii. scattering jacobi 
        # lets say it will taken through a newly introduced variable `type`.

        sysType = dict(scf.items('sysInfo'))['type']
        if sysType == 'spectroscopic': # for spectroscopic system
            self.rhoList = map(float, gInfo['rho'].split(','))
            mInfo = dict(scf.items('mInfo'))
            self.vModes = [int(i)-1 for i in mInfo['varying'].split(',')]
            if len(self.rhoList) == 1:
                raise Exception('Give a grid rho value for spectroscopic calculation')
            self.rhoGrid = self.makeGrid(self.rhoList)



        elif sysType == 'scattering_hyper': # for scattering system
            self.rhoList = map(float, gInfo['rho'].split(','))
            if len(self.rhoList) != 1:
                raise Exception('Give a fixed rho value for scattering calculation')
            self.rho = self.rhoList[0]
            self.thetaList = map(float, gInfo['theta'].split(','))
            self.thetaGrid = self.makeGrid(self.thetaList)
            if not 'scale' in self.eInfo.keys():
                raise Exception('Assymptotic energy value is mandatory for scaling the Energy surafaces for scattering system.')
            #scalling factor for scattering system

        elif sysType == 'scattering_jacobi':# for the jacobi case
            self.smallR = float(gInfo['small_r'])
            self.capR   = float(gInfo['capital_r'])
            self.gamma  = float(gInfo['gamma'])
            self.q      = float(gInfo['q'])
            # phigrid is already above

        else :
            raise Exception('Not a proper system type')


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
                gny = self.interpolate(gx,gy,grid2)
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
            m = '{:.<90}'.format(datetime.now().strftime("[%d-%m-%Y %I:%M:%S %p]     ") + m)
        else:
            m+='\n'
        self.logFile.write(m)

    def moveFiles(self, path):
        ''' Saves the geometry out and results file in a specific directory'''
        os.makedirs(path)
        files = glob('*.xyz') + glob('*.out') + glob('*.res')
        for file in files:
            shutil.move(file, path)


    def runThisMolpro(self, grid1, gridn1, grid2, gridn2, filEe, fileN1, fileN2):
        '''Runs the molpro for each gridpoints '''
        # subprocess calls blocks system I/O buffer, so the I/Os (sometimes) have to be flushed out explicitely 
        # open files to store result
        filee  = open(filEe, 'w', buffering=1)
        filen1 = open(fileN1,'w', buffering=1)
        filen2 = open(fileN2,'w', buffering=1)


        #grid1 is theta or rho
        #grid2 is phi
        #create folder for incomplete jobs , delete if already exists
        if os.path.isdir('IncompleteJobs'):
            shutil.rmtree('IncompleteJobs')
        os.mkdir('IncompleteJobs')

        if os.path.isdir('CompleteJobs'):
            shutil.rmtree('CompleteJobs')
        os.mkdir('CompleteJobs')

        self.equiRun()
        if self.__class__.__name__ == 'Spectroscopic':
            self.moveFiles('CompleteJobs/Equilibrium_point')
        else:
            self.moveFiles('CompleteJobs/Initial_point')

        done = False 

        for g2 in grid2:
            shutil.copy('molpro_init.wfu', 'molpro.wfu')      # copy the wfu from the initial/equilibrium job
            for g1 in grid1: 

                #now if rho 0 is provided run the rho=0=phi step only once
                if (g1==0) & (g2!=0) & (not done):
                    continue # don't run the job for other rhos
                else:
                    done = True

                self.createGridGeom(g1, g2)
                self.msg('Running molpro job for {} = {}, {} = {}.......'.format(gridn1, g1, gridn2, g2))
                path = 'CompleteJobs/{}_{}/{}_{}'.format(gridn1, g1, gridn2, g2)
                shutil.copy('molpro.wfu',self.scrdir)
                exitcode = subprocess.call(
                    ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'grid.com','--no-xml-output']
                    )
                if exitcode==0:
                    self.msg( ' Job successful.', cont=True)
                else:
                    self.msg(' Job failed.', cont=True)
                    path = path.replace('C', "Inc")
                    self.moveFiles(path)
                    continue

                enrData = self.parseResult('enr.res').flatten()
                tau1, tau2 = self.getTau(g1, g2)
                np.savetxt(filee,  np.append([g1,g2],enrData)[None], fmt='%.8f', delimiter='\t')
                np.savetxt(filen1, np.append([g1,g2],tau1)[None],  fmt='%.8f', delimiter='\t')
                np.savetxt(filen2, np.append([g1,g2],tau2)[None],  fmt='%.8f', delimiter='\t')
                self.moveFiles(path)
            filee.write('\n')
            filen1.write('\n')
            filen2.write('\n')
        # self.removeFiles(allOut=True)   # removes the wfu and .com files after complete run
        self.msg('All molpro jobs done.')

        scat = self.__class__.__name__=='Scattering'
        #fill the missing values by 1D interpolation
        # '_mod' suffix means files with filled data
        for file in [filEe, fileN1, fileN2]:
            dat = self.interp(file)
            dat[:,1] = np.deg2rad(dat[:,1])
            if scat:                             # for scattering 
                dat[:,0] = np.deg2rad(dat[:,0])  # for scattering also transform the column 0 i.e. theta column
            if file == filEe:                    # Now scale the energy
                if scat:                         # for scattering scale it with user given assymptote value
                    dat[:,2:] -= float(self.eInfo['scale'])
                else:                            # for spectroscopic scale with equilibrium energy value
                    dat[:,2:] -= np.loadtxt('equienr.dat')[0]
            self.writeFile(file.replace('.dat', '_mod.dat'), dat)



class Spectroscopic(Base):
    ''' Inherited from the Base class, this class contains necessary methods
     for running molpro for a Spectroscopic system'''
    def __init__(self, conFig ,atomFile, geomFile , freqFile , wilsonFile ):
        self.parseData(atomFile)
        self.parseSData(geomFile, freqFile, wilsonFile)
        self.parseConfig(conFig)

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

        curGeom  = self.equiGeom + np.einsum('ijk,k->ij',self.wilFM, qCord)
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
        exitcode = subprocess.call(
            ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'init.com','--no-xml-output']
            )
        if exitcode==0: 
            self.msg( ' Job successful', cont=True)
        else:
            self.msg( ' Job failed', cont=True)
            sys.exit('Molpro failed in equilibrium step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt='%.8f')




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
    def __init__(self, conFig, atomFile ):
        self.parseData(atomFile)
        self.parseConfig(conFig)

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
       #! do this in more short way?
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
        p1 = [0, rc*np.sin(gamma), rc*np.cos(gamma)]
        p2 = [0,0.0, -rs/2.0]
        p3 = [0,0.0 ,rs/2.0 ]
        return np.array([p1, p2, p3])*0.529177209 # return in angstrom

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
        exitcode = subprocess.call(
            ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'init.com','--no-xml-output']
            )
        if exitcode==0: 
            self.msg( ' Job successful', cont= True)
        else:
            self.msg( ' Job failed', cont= True)
            sys.exit('Molpro failed in initital step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt='%.8f')



    def runMolpro(self):
        self.runThisMolpro(
            self.thetaGrid, 
            'Theta', 
            self.phiGrid, 
            'Phi', 
            'energy.dat', 
            'tau_theta.dat', 
            'tau_phi.dat' )




class Jacobi(Base):
    def __init__(self, conFig, atomFile ):
        self.parseData(atomFile)
        self.parseConfig(conFig)


    def createDdrTemplate(self):
        ''''Creates the molpro template files ddr job'''

        cas = self.eInfo['cas']

        #remove core from the cas
        try:
            core = re.search( '[0-9a-zA-Z,;](core,\d+;*)', cas).group(1)
            cas = cas.replace(core, '').strip(';')
        except:
            pass


        #mrci has to be done for ddr for using the ddr wf in ddr nact calculation
        energyLine = """
            {{mcscf;{cas1}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;{extra2}}}
            """.format(state   =self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'])
            

        molproTemplate=textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job for analytical job.
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom1.xyz

            {enrl}

            show, energy
            table, energy
            save,enr.res,new
            {{table,____; noprint,heading}}

            !for +dr
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas1} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;{extra2}}}
            {{ci;trans,6000.2,6001.2;dm,8001.2;{extra3}}}


            !for -dr
            symmetry,nosym
            geometry=geom3.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;{extra2}}}
            {{ci;trans,6000.2,6002.2;dm,8002.2;{extra3}}}



            '''.format( memory = self.memory,
                        basis  = self.eInfo['basis'],
                        enrl   = energyLine,
                        state  = self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'],
                        extra3 = self.nInfo['nact_extra']))


        nactTemp= ''

        for i,j in self.nactPairs:
            nactTemp+=textwrap.dedent(''' 
                !for taur     
                {{ddr,{dt},2140.2,2241.2,8001.2;state, {j}.1,{i}.1}}
                nacmepv=nacme

                {{ddr,-{dt},2140.2,2242.2,8002.2;state, {j}.1,{i}.1}}
                nacmemv=nacme

                nacmp = 0.5*(nacmepv+ nacmemv)

                table, nacmp
                save,ddrnact{i}{j}.res,new;

                '''.format(dt=self.dt,dp=self.dp,i=i,j=j))



        molproTemplate += nactTemp + '\n---\n'

        with open('grid.com','w') as f:
            f.write(molproTemplate)

        molproInitTemplate = textwrap.dedent('''
            ***, Molpro template created from ADT program for ddr job.
            memory,{memory}
            file,2,molpro_init.wfu,new;

            basis={basis};

            symmetry,nosym

            geomtyp=xyz
            geometry=geom.xyz

            {{hf}}
            {enrl}

            show, energy
            table, energy
            save,equienr.res,new
            {{table,____; noprint,heading}}

            '''.format( memory = self.memory,
                        basis = self.eInfo['basis'],
                        enrl=energyLine,
                        hf = self.eInfo['hf']))

        with open('init.com', 'w') as f:
            f.write(molproInitTemplate)




    def getTauAna(self, phi):
        '''Used in Analytical NACT calculation'''
        tauph = []
        for i,j in self.nactPairs:
            file   = 'ananac{}{}.res'.format(i,j)
            grads  = self.parseResult(file)
            val = self.q*(-grads[2,0]*self.sin(phi) + grads[2,1]*self.cos(phi))
            tauph.append(np.abs(val))
        return tauph


    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        tauph = []
        for i,j in self.nactPairs : 
            file = 'ddrnact{}{}.res'.format(i,j)
            val  = self.parseResult(file)
            tauph.append(np.abs(val))
        return np.array(tauph)




    def createOneGeom(self, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given r, R, gamma, q, phi'''

        curGeom  = self.bohrtoang*np.array([
            [-self.smallR/2.0,0.0,0.0],
            [self.smallR/2.0,0.0,0.0],
            [self.capR*self.cos(self.gamma)+self.q*self.cos(phi), 
            self.capR*self.sin(self.gamma)+self.q*self.sin(phi),
            0.0]])
        msg = 'for r = {}, R = {}, gamma = {}, Phi = {}'.format(self.smallR, self.capR, self.gamma, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])
        with open(outFile,"w") as f:
            f.write(tmp)


    def createAllGeom(self, theta, phi):
        ''' Creates 3 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom(phi,    'geom1.xyz')
        self.createOneGeom(phi+self.dp,'geom2.xyz')
        self.createOneGeom(phi-self.dp,'geom3.xyz')



    def equiRun(self):
        '''Runs molpro for a initial geometry, i.e. phi=0'''
        self.createOneGeom( 0)
        self.msg( "Running molpro job for initial point....." )
        sys.stdout.flush()
        exitcode = subprocess.call(
            ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'init.com','--no-xml-output']
            )
        if exitcode==0: 
            self.msg( ' Job successful', cont= True)
        else:
            self.msg( ' Job failed', cont= True)
            sys.exit('Molpro failed in initital step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt='%.8f')



    def runMolpro(self):
        '''Runs the molpro for each gridpoints '''
        # subprocess calls blocks system I/O buffer, so the I/Os (sometimes) have to be flushed out explicitely 
        # open files to store result
        filee  = open('energy.dat', 'w', buffering=1)
        filen = open('tau_phi.dat','w', buffering=1)


        #grid1 is theta or rho
        #grid2 is phi
        #create folder for incomplete jobs , delete if already exists
        if os.path.isdir('IncompleteJobs'):
            shutil.rmtree('IncompleteJobs')
        os.mkdir('IncompleteJobs')

        if os.path.isdir('CompleteJobs'):
            shutil.rmtree('CompleteJobs')
        os.mkdir('CompleteJobs')


        self.equiRun()

        shutil.copy('molpro_init.wfu', 'molpro.wfu')
        for phi in self.phiGrid:
                self.createGridGeom(phi)
                self.msg('Running molpro job for Phi = {}.......'.format(phi))
                path = 'CompleteJobs/Phi_{}'.format(phi)
                shutil.copy('molpro.wfu',self.scrdir)

                exitcode = subprocess.call(
                    ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'grid.com','--no-xml-output']
                    )
                if exitcode==0:
                    self.msg( ' Job successful.', cont=True)

                else:
                    self.msg(' Job failed.', cont=True)
                    path = path.replace('C', "Inc")
                    self.moveFiles(path)
                    continue 
                firstJobDone = True
                enrData = self.parseResult('enr.res').flatten()
                tau = self.getTau(phi)
                self.moveFiles(path)
                np.savetxt(filee,  np.append([phi],enrData)[None], fmt='%.8f', delimiter='\t')
                np.savetxt(filen, np.append([phi],tau)[None],  fmt='%.8f', delimiter='\t')
        self.msg('All molpro jobs done.', cont=True)

        for file in ['energy.dat', 'tau_phi.dat']:
            dat = np.loadtxt(file)
            grid = self.phiGrid
            res = grid
            for col in range(1, dat.shape[1]):
                filledDat = self.interpolate(dat[:,0], dat[:,col], grid)
                res = np.column_stack([res, filledDat]) 
            res[:,0] = np.deg2rad(res[:,0]) # convert phi to radian
            if file == 'energy.dat':
                res[:,1:] -= float(self.eInfo['scale'])
            np.savetxt(
                file.replace('.dat', '_mod.dat'),
                res, fmt='%.8f', delimiter='\t'
            )




# if __name__ == "__main__":
#     s = Jacobi('./molpro.config', './atomfile.dat')
#     s.createOneGeom(0)
#     #s = Spectroscopic('./molpro.config', './atomfile.dat', 'geomfile.dat','frequency.dat', 'wilson.dat')
#     s.runMolpro()
