from __future__ import absolute_import, unicode_literals, division, print_function
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
from glob import glob
from datetime import datetime
from adt.numeric.adtmod import adt as fadt


if sys.version_info.major>2:
    from configparser import ConfigParser as ConfigParser
else :
    from ConfigParser import SafeConfigParser as ConfigParser

# 180/pi multiplication from ddr removed


def mainFunction(logger, conFig, atomFile, *args):
    '''Read the configuration file and run corresponding calculation'''
    scf = ConfigParser()
    scf.read(conFig)
    sysType = scf.get('sysInfo', 'type')

    txt = "Starting molpro jobs."

    if sysType == 'spec':
        # args is a list of [geomFile, freqfile, wilsonFile]
        jobRunner = Spectroscopic(scf, atomFile, *args)
        txt+='''
            System type               : Spectroscopic
            Co-ordinate type          : Normal Modes
            Molpro Config file        : {}
            Atom Info file            : {}
            Equilibrium Geometry file : {}
            Frequency Info file       : {}
            Wilson Matrix file        : {}'''.format(conFig, atomFile, *args)

    elif sysType == 'scat_hyper':
        jobRunner = Scattering(scf, atomFile)
        txt+= '''
            System type               : Scattering
            Co-ordinate type          : Hyperspherical
            Molpro Config file        : {}
            Atom Info file            : {}'''.format(conFig, atomFile)

    elif sysType == 'scat_jacobi':
        jobRunner = Jacobi(scf, atomFile)
        txt+='''
            System type               : Scattering
            Co-ordinate type          : Jacobi
            Molpro Config file        : {}
            Atom Info file            : {}'''.format(conFig, atomFile)

    else :
        raise Exception('Not a proper system type')
    logger.info(txt+"\nCheck 'adt_molpro.log' for progress.\n")
    fls = jobRunner.runMolpro()
    txt = """Molpro jobs completed. \n\tData saved in the following files
        Energy file : energy_mod.dat"""
    if len(fls)==3: # 2D case
        txt+="""
        NACT 1 file : {}
        NACT 2 file : """.format(fls[1])
    else :         # 1D jacobi case
        txt+="""
        NACT file   : """
    txt+=fls[-1] +'\n\n'
    logger.info(txt)

    return fls



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
        diff = fadt.spline(x,y)
        return np.array([ fadt.splint(x,y,diff,xi) for xi in newx])

    def createAnaTemplate(self):
        ''''Creates the molpro template files analytical job'''

        # remove core from cas for multi calculation
        cas = re.sub('[0-9a-zA-Z,;](core,\d+;*)', '', self.eInfo['cas'])


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
                nactTemp += "cpmcscf,nacm,{f}.1,{s}.1,record=51{n:02}.1,{extra};\n".format(
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


        cas = re.sub('[0-9a-zA-Z,;](core,\d+;*)', '', self.eInfo['cas'])

        #mrci has to be done for ddr nact calculation
        energyLine = """
            {{mcscf;{cas1}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;dm,8000.2;{extra2}}}
            """.format(state   =self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'])

        # d1 is increment in rho, theta or q depending on calculation type
        # and d2 is always phi
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

            !for +d1
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas1} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;{extra2}}}
            {{ci;trans,6000.2,6001.2;dm,8001.2;{extra3}}}


            !for -d1
            symmetry,nosym
            geometry=geom3.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2242.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6002.2;{extra2}}}
            {{ci;trans,6000.2,6002.2;dm,8002.2;{extra3}}}



            !for +d2
            symmetry,nosym
            geometry=geom4.xyz
            {{multi;{cas1}; wf,{wf};state,{state};start,2140.2;orbital,2243.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6003.2;{extra2}}}
            {{ci;trans,6000.2,6003.2;dm,8003.2;{extra3}}}


            !for -d2
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

            #implementing three point cebtral difference
            nactTemp+=textwrap.dedent(''' 
                !for taur     
                {{ddr, 2*{d1}
                orbital,2140.2,2141.2,2142.2;
                density,8000.2,8001.2,8002.2;
                state, {j}.1,{i}.1
                }}
                nacmr = nacme

                !for taup
                {{ddr, 2*{d2}
                orbital,2140.2,2143.2,2144.2;
                density,8000.2,8003.2,8004.2;
                state, {j}.1,{i}.1
                }}
                nacmp = nacme

                table, nacmr,nacmp
                save,ddrnact{i}{j}.res,new;

                '''.format(d1=self.d1,d2=self.d2,i=i,j=j))


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
        # due to precision error from np.arange, linspace is used
        bl = int(round((ll[1]-ll[0])/ll[2]+1))
        return np.linspace(ll[0], ll[1], bl)

    def parseData(self, atomFile):
        ''' Parse atom names and masses from data file'''
        atomData = np.loadtxt(atomFile, 
            dtype={'names': ('names', 'mass'),'formats': ('U1', 'f8')})
        self.atomNames = atomData['names']
        self.atomMass  = atomData['mass']

    def parseCommonConfig(self, scf):
        ''' 
        Parses configuration for running molpro from the provided molpro config file 
        and sets up different methods and attributes relavant to the configuration 
        '''

        molInfo =  dict(scf.items('molInfo'))
        self.scrdir = molInfo['scrdir']
        self.memory = molInfo['memory']
        try:
            self.proc = molInfo['processor']
        except KeyError:
            self.proc = '1'

        self.eInfo = dict(scf.items('eInfo'))
        self.nInfo = dict(scf.items('nInfo'))
        self.gInfo = dict(scf.items('gInfo'))

        #set some default values for the optional arguments
        if not 'multi_extra' in self.eInfo.keys():
            self.eInfo['multi_extra'] = ''
        if not 'mrci_extra' in self.eInfo.keys():
            self.eInfo['mrci_extra'] = ''
        if not 'nact_extra' in self.nInfo.keys():
            self.nInfo['nact_extra'] = ''


        # set up initial HF calculation according to 'restricted' parameter
        self.eInfo['hf'] = 'uhf'
        try:
            if self.eInfo['restricted'].lower() == 'true':
                self.eInfo['hf'] = 'hf'
        except:
            pass
        if 'uhf_extra' in self.eInfo.keys():
            self.eInfo['hf'] += ';' + self.eInfo['uhf_extra']


        self.state = int(self.eInfo['state'])
        self.nactPairs = [[i, j] for j in range(2, self.state+1) for i in range(1, j)]

        self.nTau = len(self.nactPairs)


        self.logFile = open('adt_molpro.log', 'w')

    def parseResult(self, file):
        ''' Parses result from the ouptpu .res files'''
        with open(file,"r") as f:
            dat = f.read().replace("D","E").strip().split("\n")[1:]
        # dat = [float(j) for i in dat for j in i.strip().split()]#[map(float,i.strip().split()) for i in dat]
        dat = [[float(j) for j in i.strip().split()] for i in dat]
        return np.array(dat)

    def writeFile(self, file, data):
        ''' Writes output data in plain txt'''
        file = open(file,'w')
        for tp in np.unique(data[:,0]):
            np.savetxt( file, data[data[:,0]==tp] ,delimiter="\t", fmt=str("%.8f"))
            file.write("\n")
        file.flush()

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
        self.logFile.flush()

    def moveFiles(self, path):
        ''' Saves the geometry out and results file in a specific directory'''
        os.makedirs(path)
        files = glob('*.xyz') + glob('*.out') + glob('*.res') + glob('*.xml')
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
                    ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'grid.com']
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
                np.savetxt(filee, np.append([g1,g2],enrData)[None], fmt=str("%.8f"), delimiter='\t')
                np.savetxt(filen1, np.append([g1,g2],tau1)[None],  fmt=str("%.8f"), delimiter='\t')
                np.savetxt(filen2, np.append([g1,g2],tau2)[None],  fmt=str("%.8f"), delimiter='\t')
                self.moveFiles(path)
            filee.write('\n')
            filen1.write('\n')
            filen2.write('\n')
        # self.removeFiles(allOut=True)   # removes the wfu and .com files after complete run
        self.msg('All molpro jobs done.\n' ,cont=True)





        scat = self.__class__.__name__=='Scattering'  # check which class is calling this
        # only for scattering case the 1st column i.e. thete need to be converted to radian
        #fill the missing values by 1D interpolation
        # '_mod' suffix means files with filled data
        inFiles = [filEe, fileN1, fileN2]
        outFiles = [i.replace('.dat', '_mod.dat') for i in inFiles]
        for iFile,oFile in zip(inFiles, outFiles):
            dat = self.interp(iFile)
            dat[:,1] = np.deg2rad(dat[:,1])      # second column is radian, convert it to phi
            if scat:                             # for scattering 
                dat[:,0] = np.deg2rad(dat[:,0])  # for scattering also transform the column 0 i.e. theta column
            if iFile == filEe:                    # Now scale the energy
                if scat:                         # for scattering scale it with user given assymptote value
                    dat[:,2:] -= float(self.eInfo['scale'])
                else:                            # for spectroscopic scale with equilibrium energy value
                    dat[:,2:] -= np.loadtxt('equienr.dat')[0]
            self.writeFile(oFile, dat)



class Spectroscopic(Base):
    ''' Inherited from the Base class, this class contains necessary methods
     for running molpro for a Spectroscopic system'''
    def __init__(self, conFig ,atomFile, geomFile , freqFile , wilsonFile ):
        self.parseData(atomFile)
        self.parseSData(geomFile, freqFile, wilsonFile)
        self.parseCommonConfig(conFig)
        self.parseThisConfig(conFig)


    def parseThisConfig(self, conFig):
        '''Parses some extra configuaration specific to this type'''

        mInfo = conFig.get('mInfo', 'varying')
        self.vModes = [int(i)-1 for i in mInfo.split(',')]

        self.rhoList = [float(i) for i in self.gInfo['rho'].split(',')]
        assert len(self.rhoList)==3, "A grid of rho is required"
        self.rhoGrid = self.makeGrid(self.rhoList)

        self.phiList = [float(i) for i in self.gInfo['phi'].split(',')]
        assert len(self.phiList)==3, "A grid of phi is required"
        self.phiGrid = self.makeGrid(self.phiList)

        if self.nInfo['method']=='cpmcscf':
            self.createAnaTemplate()
            self.getTau = self.getTauAna
            self.createGridGeom = self.createOneGeom


        elif self.nInfo['method']=='ddr':
            self.d1 = float(gInfo['drho'])
            self.d2 = float(gInfo['dphi'])
            self.createDdrTemplate()
            self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr



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
        tau =  np.stack([self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs]).T
        return np.abs(tau)



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
            ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'init.com']
            )
        if exitcode==0: 
            self.msg( ' Job successful', cont=True)
        else:
            self.msg( ' Job failed', cont=True)
            sys.exit('Molpro failed in equilibrium step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt=str("%.8f"))




    def runMolpro(self):
        file = ['energy.dat','tau_rho.dat','tau_phi.dat']
        self.runThisMolpro(
            self.rhoGrid, 
            'Rho', 
            self.phiGrid, 
            'Phi', *file )
        return [i.replace('.dat', '_mod.dat') for i in file]



class Scattering(Base):
    ''' Inherited from the Base class, this class containg necessary methods
     for running molpro for a Scattering system in hyperspherical coordinate'''
    def __init__(self, conFig, atomFile):
        self.parseData(atomFile)
        self.parseCommonConfig(conFig)
        self.parseThisConfig(conFig)


    def parseThisConfig(self, conFig):
        self.rhoList = [float(i) for i in self.gInfo['rho'].split(',')]
        assert len(self.rhoList) ==1 , 'A fixed rho value is required'
        self.rho = self.rhoList[0]

        self.thetaList = [float(i) for i in self.gInfo['theta'].split(',')]
        assert len(self.thetaList)==3, "A grid of phi is required"
        self.thetaGrid = self.makeGrid(self.thetaList)

        self.phiList = [float(i) for i in self.gInfo['phi'].split(',')]
        assert len(self.phiList)==3, "A grid of phi is required"
        self.phiGrid = self.makeGrid(self.phiList)

        if not 'scale' in self.eInfo.keys():
            raise Exception('Assymptotic energy value is mandatory for scaling the Energy surafaces for scattering system.')


        if self.nInfo['method']=='cpmcscf':
            self.createAnaTemplate()
            self.getTau = self.getTauAna
            self.createGridGeom = self.createOneGeom

        elif self.nInfo['method']=='ddr':
            self.d1 = float(self.gInfo['dtheta'])
            self.d2 = float(self.gInfo['dphi'])
            self.createDdrTemplate()
            self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr


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
        eps3 = np.rad2deg(eps3)
        eps2 = np.rad2deg(eps2)
        R1 = (1.0/np.sqrt(2.0))*self.rho*d3*np.sqrt(1.0+ self.sin(theta)*self.cos(phi+eps3)) 
        R2 = (1.0/np.sqrt(2.0))*self.rho*d1*np.sqrt(1.0+ self.sin(theta)*self.cos(phi))      
        R3 = (1.0/np.sqrt(2.0))*self.rho*d2*np.sqrt(1.0+ self.sin(theta)*self.cos(phi-eps2)) 

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
        self.createOneGeom(theta+self.d1, phi,   'geom2.xyz')
        self.createOneGeom(theta-self.d2, phi,   'geom3.xyz')
        self.createOneGeom(theta,    phi+self.d2,'geom4.xyz')
        self.createOneGeom(theta,    phi-self.d2,'geom5.xyz')

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
        dTheta = self.thetaList[2]/100.0
        dPhi = self.phiList[2]/100.0

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
        tau = np.vstack([self.parseResult('ddrnact{}{}.res'.format(i,j)) 
                                    for i,j in self.nactPairs]).T
        return np.abs(tau)


    def equiRun(self):
        '''Runs molpro for a initial geometry, i.e. theta=phi=0'''
        self.createOneGeom(0, 0)
        self.msg( "Running molpro job for initial point....." )
        sys.stdout.flush()
        exitcode = subprocess.call(
            ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'init.com']
            )
        if exitcode==0: 
            self.msg( ' Job successful', cont= True)
        else:
            self.msg( ' Job failed', cont= True)
            sys.exit('Molpro failed in initital step')
        equiData = self.parseResult('equienr.res').flatten()
        np.savetxt('equienr.dat', equiData, fmt=str("%.8f"))



    def runMolpro(self):
        file = ['energy.dat','tau_theta.dat','tau_phi.dat']
        self.runThisMolpro(
            self.thetaGrid, 
            'Theta', 
            self.phiGrid, 
            'Phi', *file )
        return [i.replace('.dat', '_mod.dat') for i in file]






class Jacobi(Base):
    def __init__(self, conFig, atomFile):
        self.parseData(atomFile)
        self.parseCommonConfig(conFig)
        self.parseThisConfig(conFig)


    def parseThisConfig(self, conFig):
        self.phiList = [float(i) for i in self.gInfo['phi'].split(',')]
        assert len(self.phiList)==3, "A grid of phi is required"
        self.phiGrid = self.makeGrid(self.phiList)

        self.smallR = float(self.gInfo['small_r'])
        self.capR   = float(self.gInfo['capital_r'])
        self.gamma  = float(self.gInfo['gamma'])
        # Checking if q is given in a list for grid or just
        # a single value of q depending on that 2D or 1D adt will be done
        ql = self.gInfo['q'].split(',')
        if len(ql) == 1:  # 1D Jacobi will be done
            self.Jacobi1D          = True
            self.q                 = float(self.gInfo['q'])
            self.runMolpro         = self.runMolpro1D
            self.createDdrTemplate = self.createDdrTemplate1D
            self.createOneGeom     = self.createOneGeom1D
            self.createAllGeom     = self.createAllGeom1D
        else :            # 2D Jacobi will be done
            self.Jacobi1D = False
            ql = [float(i) for i in ql]
            self.qGrid = self.makeGrid(ql)

        if not 'scale' in self.eInfo.keys():
            raise Exception('Assymptotic energy value is mandatory for scaling the Energy surafaces for scattering system.')


        if self.nInfo['method']=='cpmcscf':
            self.createAnaTemplate()
            self.getTau = self.getTauAna
            self.createGridGeom = self.createOneGeom

        elif self.nInfo['method']=='ddr':
            if Jacobi1D:
                self.d1 = float(gInfo['dphi'])
            else:
                self.d1 = float(gInfo['dq'])
                self.d2 = float(gInfo['dphi'])
            self.createDdrTemplate()
            self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr



    def createDdrTemplate1D(self):
        ''''Creates the molpro template files ddr job'''

        cas = re.sub( '[0-9a-zA-Z,;](core,\d+;*)', '', self.eInfo['cas'])


        #mrci has to be done for ddr for using the ddr wf in ddr nact calculation
        energyLine = """
            {{mcscf;{cas1}; wf,{wf};state,{state};start,2140.2; orbital,2140.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6000.2;dm,8000.2;{extra2}}}
            """.format(state   =self.eInfo['state'],
                        wf     = self.eInfo['wf'],
                        cas    = self.eInfo['cas'],
                        cas1   = cas,
                        extra1 = self.eInfo['multi_extra'],
                        extra2 = self.eInfo['mrci_extra'])
            
        # for jacobi 1d d1 is phi, unlike 2d where d1 is q.
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

            !for +d1
            symmetry,nosym
            geometry=geom2.xyz
            {{multi;{cas1} ;wf,{wf};state,{state};start,2140.2;orbital,2241.2;{extra1}}}
            {{mrci; {cas}; wf,{wf};state,{state};save,6001.2;{extra2}}}
            {{ci;trans,6000.2,6001.2;dm,8001.2;{extra3}}}


            !for -d1
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
                {{ddr, 2*{d1}
                orbital,2140.2,2141.2,2142.2;
                density,8000.2,8001.2,8002.2;
                state, {j}.1,{i}.1
                }}

                table, nacme
                save,ddrnact{i}{j}.res,new;

                '''.format(d1=self.d1,i=i,j=j))


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
        return (180.0/np.pi)*np.array(tauph)


    def createOneGeom(self, q, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given r, R, gamma, q, phi'''

        curGeom  = self.bohrtoang*np.array([
            [-self.smallR/2.0,0.0,0.0],
            [self.smallR/2.0,0.0,0.0],
            [self.capR*self.cos(self.gamma)+q*self.cos(phi), 
            self.capR*self.sin(self.gamma)+q*self.sin(phi),
            0.0]])
        msg = 'for r = {}, R = {}, gamma = {}, Phi = {}, q= {}'.format(self.smallR, self.capR, self.gamma, phi, q)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])
        with open(outFile,"w") as f:
            f.write(tmp)


    def createAllGeom(self, q, phi):
        ''' Creates 5 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom(q, phi, 'geom1.xyz')
        self.createOneGeom(q+self.q, phi, 'geom2.xyz')
        self.createOneGeom(q-self.q, phi, 'geom3.xyz')
        self.createOneGeom(q, phi+self.dp,'geom4.xyz')
        self.createOneGeom(q, phi-self.dp,'geom5.xyz')



    def createOneGeom1D(self, phi, outFile='geom.xyz'):
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


    def createAllGeom1D(self, phi):
        ''' Creates 3 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom1D(phi,    'geom1.xyz')
        self.createOneGeom1D(phi+self.dp,'geom2.xyz')
        self.createOneGeom1D(phi-self.dp,'geom3.xyz')


    def runMolpro(self):
        file = ['energy.dat','tau_q.dat','tau_phi.dat']
        self.runThisMolpro(
            self.thetaGrid, 
            'q', 
            self.phiGrid, 
            'Phi', *file )
        return [i.replace('.dat', '_mod.dat') for i in file]



    def runMolpro1D(self):
        '''Runs the molpro for each gridpoints '''
        # this here is only for 1d jacobi
        # subprocess calls blocks system I/O buffer, so the I/Os (sometimes) have to be flushed out explicitely 
        # open files to store result
        filee  = open('energy.dat', 'w', buffering=1)
        filen  = open('tau_phi.dat','w', buffering=1)


        if os.path.isdir('IncompleteJobs'):
            shutil.rmtree('IncompleteJobs')
        os.mkdir('IncompleteJobs')

        if os.path.isdir('CompleteJobs'):
            shutil.rmtree('CompleteJobs')
        os.mkdir('CompleteJobs')
        firstDone = False

        # shutil.copy('molpro_init.wfu', 'molpro.wfu')
        for phi in self.phiGrid:
                self.createGridGeom(phi)
                self.msg('Running molpro job for Phi = {}.......'.format(phi))
                path = 'CompleteJobs/Phi_{}'.format(phi)
                if firstDone:
                    shutil.copy('molpro.wfu',self.scrdir)

                exitcode = subprocess.call(
                    ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'grid.com']
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
                np.savetxt(filee,  np.append([phi],enrData)[None], fmt=str("%.8f"), delimiter='\t')
                np.savetxt(filen, np.append([phi],tau)[None],  fmt=str("%.8f"), delimiter='\t')
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
                res, fmt=str("%.8f"), delimiter='\t'
            )
        return ['energy_mod.dat', 'tau_phi_mod.dat']

