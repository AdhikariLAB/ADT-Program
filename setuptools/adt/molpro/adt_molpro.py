
# wf now splitted into electron, spin and charge numbers and the wf card will be built inside the scripts
#? the symmetry information is inside the sysInfo section ????
# so eInfo now has this mandatory keyworsd

# 1. method    ----> same as before
# 2. basis     ----> same as before
# 3. cas       ----> occ,2;closed,2         no symmetry
#                    occ,2,3;closed,2,3     Cs
#                    occ,2,3,4,5;...........C2v so on...
# 4. electron  ----> number of electron, one number
# 5. spin      ----> one number
# 6. charge    ----> one number
# 7. state     ----> n numbers seperated by comma, n= number of IREPS 
#                    0 is mandatory for no state calculation for a particular state



###### An these optional keywords
# 1. scale 
# 2. uhf extra     ---> what will be done with this when wf is provided through this ?????
# 3. multi extra 
# 4. mrci extra






#* all the runMolpro function, after running the  ab initio jobs returns two things 
#* (i) a string indicating 1D or 2D ab initio and 
#* (ii) a list contiainig files  that will be used in ADT calculation in list of different IREPs

# ! IREP for the wf card for uhf extra???







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
    typ, files = fls
    txt = """Molpro jobs completed. Results stored in following files"""


    files = [', '.join(i) for i in zip(*filter(None, files))]
    if typ=="2D":
        txt+= """
            Energy files : {},
            NACT 1 files : {},
            NACT 2 files : {}
        """.format(*files)
    else: # 1D jacobi case
        txt+= """
            Energy files : {},
            NACT 1 files : {}
        """.format(*files)

    logger.info(txt)

    return fls








class Base(object):
    '''
        A base class containing the common methods to be used in both cases of spectroscopic and scattering
    '''
    # some constants to be used later
    angtobohr    = 1.8897259886
    hbar         = 0.063508
    cInvToTauInv = 0.001883651
    bohrtoang    = 0.529177

    # this dictionary keeps track of number of IREP for a particualr symmetry
    # want to include new symmetry, just update this.
    IREP_LIST = {
        "x": 2,
        "nosym" : 1
    }


    def sin(self, x):
        # A sin function that directly takes degree as input unlike numpy
        return np.sin(np.deg2rad(x))

    def cos(self, x):
        # A cos function that directly takes degree as input unlike numpy
        return np.cos(np.deg2rad(x))

    def interpolate(self, x, y, newx):
        # an function to calculate 1d interpolation of the new grid `newx` from older `x,y`
        # This is imported from the fortran 1D interpolation subroutine from numeric section as
        # numpy does'nt have any cubic spline ineterpolation.
        diff = fadt.spline(x,y)
        return np.array([ fadt.splint(x,y,diff,xi) for xi in newx])



    def irepChecks(self):
        # Checks everythings against the IREP numbers

        
        self.nStateList = [int(i) for i in re.findall('(\d+)', self.eInfo['state'])]

        #* sanity checks for state
        assert len(self.nStateList)==self.nIREPs, (
                "{a} number of states required for {a} IREPs of {c} symmetry".format(
                                            a = self.nIREPs,
                                            c = self.symmetry
                                            )
                )



        # nactPairsList is a  3 dimensional list. i.e list containing the list of pairs
        self.nactPairsList =[[[i, j] for j in range(2, state+1) for i in range(1, j)] for state in self.nStateList]
        self.nTau = [len(item) for item in self.nactPairsList]


        # Chekc occ closed and core for against IREP number
        for substring in self.eInfo['cas'].split(';'):
            for item in ['occ', 'closed', 'core']:
                if item in substring:
                    assert len(re.findall('(\d+)', substring)) == self.nIREPs, (
                        "{a} number of {b} required for {a} IREPs of {c} symmetry".format(
                                                    a = self.nIREPs,
                                                    b = item,
                                                    c = self.symmetry
                                                    )
                        )


    def getEnergy(self):
        eList = np.array([])
        if (self.eInfo['method'] == 'mrci') or (self.nInfo['method']=='ddr') :
            for nIrep, state in enumerate(self.nStateList, start=1):
                if(state):
                    file = 'enr.{}.res'.format(nIrep)
                    res = self.parseResult(file).flatten()
                    eList =np.append(eList, res)
        else:
            res = self.parseResult('enr.res').flatten()
            eList =np.append(eList, res)
        return eList





    def createTemplate(self):
        # Create three molpro template -> 1. init.com, 2. grid.com and 3. scale.com, used for running molpro

        # check everything is on par with the number of irep
        self.irepChecks()

        # remove `core` from cas for mcscf
        cas = re.sub('[0-9a-zA-Z,;](core[0-9,]+;?)','', self.eInfo['cas'])

        irepCards = []
        # crate the ptoper wf keyword 
        for nIrep, state in enumerate(self.nStateList, start=1):
            if(state):
                irepCards.append( 'wf,{elec},{irep},{spin},{chrg};state,{stat}'.format(
                                elec = self.eInfo['electron'],
                                irep = nIrep,
                                spin = self.eInfo['spin'],
                                chrg = self.eInfo['charge'],
                                stat = state
                ))


        energyLine =textwrap.dedent('{{mcscf;{cas};{irepCard};start,2140.2; orbital,2140.2;{extra}}}\n'.format(
                            cas = cas,   # <<<----- the cas after stripping the `core` from the cas provided
                            irepCard =  ';'.join(irepCards),
                            extra= self.eInfo['multi_extra']
                            ))

        # If mrci is provided or ddr is to be done then do the mrci line


        # if basis is changed before nact during ddr then another mrci has to be done so fon't save the dm
        if ((self.eInfo['method'] == 'mrci') and (self.nInfo['method']!='ddr')) or((self.nInfo['method']=='ddr') and (self.nInfo['basis'] !=self.eInfo['basis'])) : #mrci but no ddr
            for irepCard in irepCards:
                irep = irepCard.split(',')[2]
                energyLine+=textwrap.dedent('''
                            {{mrci;{cas};{irepCard};{extra}}}

                            show, energy
                            table, energy
                            save,enr.{irep}.res,new
                            {{table,____; noprint,heading}}
                            
                            '''.format(
                                    cas=self.eInfo['cas'],
                                    irepCard = irepCard,
                                    extra = self.eInfo['mrci_extra'],
                                    irep = irep
                                    ))

        elif (self.nInfo['method']=='ddr') and (self.nInfo['basis'] ==self.eInfo['basis']) :
            for irepCard in irepCards:
                irep = int(irepCard.split(',')[2]) # get the irep number from the irep card
                # the saved orbital and dm depend on irep
                energyLine+=textwrap.dedent('''
                            {{mrci;{cas};{irepCard};{extra};save,60{irep:02}.2;dm,80{irep:02}.2;}}

                            show, energy
                            table, energy
                            save,enr.{irep}.res,new
                            {{table,____; noprint,heading}}
                            
                            '''.format(
                                    cas=self.eInfo['cas'],
                                    irepCard = irepCard,
                                    extra = self.eInfo['mrci_extra'],
                                    irep = irep
                                    ))

        else : # no ddr, no mrci
            #! CAUTION: for this case structure of energy file is different, remember that while reading
            energyLine += textwrap.dedent('''
                        show, energy
                        table, energy
                        save,enr.res,new
                        {table,____; noprint,heading}
            
            ''')



        # proper template will be created from this, after removing the relavant parts
        baseTemplate = textwrap.dedent('''
            ***, Molpro template created from ADT program for analytical job
            memory,{memory}
            file,2,molpro.wfu;

            basis={basis};

            symmetry,{sym}

            geometry=geom.xyz
            !HF:{{hf}}


            '''.format(memory = self.memory,
                       basis  = self.eInfo['basis'],
                       sym    = self.symmetry,
                       hf     = self.eInfo['hf']))

        baseTemplate += energyLine

        # remove the proper wildcards for proper template 
        scaleTemplate= baseTemplate.replace('molpro.wfu', 'molpro_init.wfu,new') + '\n---' # template for scaling job
        initTemplate = scaleTemplate.replace('!HF:', '')   # open the hf card for initial job

        gridTemplate = re.sub("!HF:.*", '', baseTemplate) # remove hf card, not actually necessary, as its already commented




        # Crete nact calculation part of the template 
        if self.nInfo['method']=='ddr' :
            gridTemplate += self.ddrTemplate(irepCards)
        elif self.nInfo['method']=='cpmcscf':
            gridTemplate += self.anaTemplate(irepCards)
        else:
            raise Exception("%s is not a valid NACT calculation method"%self.nInfo['method'])


        # write the templates
        with open('init.com', 'w') as f: f.write(initTemplate)

        with open('scale.com','w') as f: f.write(scaleTemplate)

        with open('grid.com', 'w') as f: f.write(gridTemplate)



    def anaTemplate(self, irepCards):
        # Create NACT part of the template for mcscf nact calculation
        # remove `core` from cas for mcscf
        cas = re.sub('[0-9a-zA-Z,;](core[0-9,]+;?)','', self.eInfo['cas'])

        # returns a list in chunks of 5 of the input sequence
        getChunk = lambda seq: [seq[i:i+5] for i in range(0, len(seq), 5)]


        # Here's the tricky part all mcscf has to be done in single line including all of the irep
        # but cpmcscf can handle only 5 at a time, so create a list of irep and pairs and loop through
        # 5 chunks at a time
        cpmcList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for f,s in nactPairs:
                cpmcList.append([ nIrep, f, s])


        fullIrepCard = ';'.join(irepCards)


        nactText = "\n\nbasis={}\n".format(self.nInfo['basis'])

        for pairChunk in getChunk(cpmcList):
            # for this one chunk of 5 pairs, mcscf has to be done once, including all the ireps
            # But those cpmcscf has to be done seperately, taking care of their ireps
            nactText += "\n{{mcscf;{cas}; {irepCard}; start,2140.2;{extra};\n".format(
                    cas=cas,
                    irepCard = fullIrepCard,
                    extra= self.eInfo['multi_extra'])


            forceText = ''
            for count, (nIrep, f,s) in enumerate(pairChunk):
                nactText += "cpmcscf,nacm,{f}.{irep},{s}.{irep},record=51{n:02}.1,{extra};\n".format(
                            f=f, s=s, 
                            n=count, 
                            irep = nIrep,
                            extra=self.nInfo['nact_extra'])

                forceText +=textwrap.dedent("""
                            force;nacm,51{n:02}.1;varsav
                            table,gradx,grady,gradz;
                            save,ananac_{a}.{r}_{b}.{r}.res,new;
                            {{table,____;  noprint,heading}}

                            """.format(a=f, b=s,r=nIrep, n=count))
            nactText += "}\n\n" + forceText

        return nactText + '\n---\n'




    def ddrTemplate(self, irepCards):
        # Create NACT part of the template for ddr nact calculation
        # d1 is increment in rho, theta or q depending on calculation type
        # and d2 is always phi
        # remove `core` from cas for mcscf
        cas = re.sub('[0-9a-zA-Z,;](core[0-9,]+;?)','', self.eInfo['cas'])

        # NOTE: a wild card `!N1D:` for the 1D ddr template purpose. 
        # When to create ddr nact template for 1D contour ab initio, 
        # we will just remove blocks following the wildcard


        # to be used inside mcscf
        fullIrepCard = ';'.join(irepCards)
        nactText = ''
        mcOrb = '2140.2'

        # if basis is changed for ddr then, do all the mcscf and mrci for the origianl geometry
        if (self.nInfo['basis'] !=self.eInfo['basis']):
            mcOrb = '2150.2'
            nactText = 'basis = {}\n'.format(self.nInfo['basis'])
            nactText+='{{mcscf;{cas};{irepCard};start,2140.2; orbital,2150.2;{extra}}}\n'.format(
                            cas = cas,   # <<<----- the cas after stripping the `core` from the cas provided
                            irepCard =  fullIrepCard,
                            extra= self.eInfo['multi_extra']
                            )
            for irepCard in irepCards:
                # the saved orbital and dm depend on irep
                irep = int(irepCard.split(',')[2])
                nactText+='{{mrci;{cas};{irepCard};{extra};save,60{irep:02}.2;dm,80{irep:02}.2;}}\n'.format(
                                    cas=self.eInfo['cas'],
                                    irepCard = irepCard,
                                    extra = self.eInfo['mrci_extra'],
                                    irep = irep
                                    )


        # How this works?
        # The index loop, loops through the cases +d1, -d1, +d2, -d2, associated with index value 0,1,2 and 3.
        # for a specific index, load appropriate geometry, do mcscf once with all ireps and save at orbital 2241,2242,2243,2244
        # for a specific index/geometry do mrci, in loop of irep and save the wavefunciotn and dm
        # for +d1 save -> 6101 onwards, dm -> 8101 onwards
        # for -d1 save -> 6126 onwards, dm -> 8126 onwards
        # for +d2 save -> 6151 onwards, dm -> 8151 onwards
        # for -d2 save -> 6176 onwards, dm -> 8176 onwards
        # so clearly this approach doesn't work if anyone is doing ddr with more than 25 irep, thats a more than safe assumption
        # the `!N1D` wild card is there to assist it in creating 1D ddr template

        indexText = {
                        0: '\n\n!for +d1',
                        1: '\n\n!for -d1',
                        2: '\n\n!N1DB\n!for +d2',
                        3: '\n\n!for -d2'
                    }


        for index in range(4):
            nactText += indexText[index]

            nactText+=textwrap.dedent(''' 
                        symmetry,{sym}
                        geometry=geom{ge:1}.xyz
                        {{multi;{cas};{irepCard};start,{mcOrb};orbital,214{ml:1}.2;{extra}}}
            '''.format(
                        irepCard = fullIrepCard,
                        cas   = cas,
                        mcOrb = mcOrb,
                        ml = index+1,
                        ge = index+2,
                        sym    = self.symmetry,
                        extra = self.eInfo['multi_extra'],
                    ))

            for irepCard in irepCards:
                irep = int(irepCard.split(',')[2])
                nactText +=textwrap.dedent(''' 
                        {{mrci;{cas};{irepCard};save,61{irep2:02}.2;{extra2}}}
                        {{ci;trans,60{irep1:02}.2,61{irep2:02}.2;dm,81{irep2:02}.2;{extra3}}}
                '''.format(
                    cas = self.eInfo['cas'],
                    irepCard = irepCard,
                    extra2 = self.eInfo['mrci_extra'],
                    extra3 = self.nInfo['nact_extra'],
                    irep1 = irep,
                    irep2 = irep + index*25
                ))

        # an wild card
        nactText +="!N1DB\n"

        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for i,j in nactPairs:
                #implementing three point cebtral difference for ddr nact calculation
                nactText+=textwrap.dedent(''' 
                    !for taur     
                    {{ddr, 2*{d1}
                    orbital,{mcOrb},2141.2,2142.2;
                    density,80{irep:02}.2,81{irep1:02}.2,81{irep2:02}.2;
                    state, {j}.{irep},{i}.{irep}
                    }}
                    nacmr = nacme


                    !N1D:for taup
                    {{ddr, 2*{d2}
                    orbital,{mcOrb},2143.2,2144.2;
                    density,80{irep:02}.2,81{irep3:02}.2,81{irep4:02}.2;
                    state, {j}.{irep},{i}.{irep}
                    }}
                    nacmp = nacme

                    table, nacmr,nacmp
                    save,ddrnact_{i}.{irep}_{j}.{irep}.res,new;

                    '''.format(
                        d1=self.d1,
                        d2=self.d2,
                        i=i,j=j,
                        mcOrb = mcOrb,
                        irep=nIrep,
                        irep1 = nIrep,
                        irep2 = nIrep+25,
                        irep3 = nIrep+50,
                        irep4 = nIrep+75,
                        ))


        return nactText + '\n---\n'



    def makeGrid(self, ll):
        """
            Create grid from a given list
        """
        # ll is the list of [start, end, step]
        # easiest way may be to use arange of numpy, 
        # but due to precison error, arange sometimes cant handle this properly and
        # sometimes skip a value the following is just a workaround to prevent that, 
        # don't know if this can circumvent that problem all the time
        bl = int(round((ll[1]-ll[0])/ll[2]+1))
        return np.linspace(ll[0], ll[1], bl)


    def parseData(self, atomFile):
        ''' 
            Parse atom names and masses from data file
        '''
        atomData = np.loadtxt(atomFile, 
            dtype={'names': ('names', 'mass'),'formats': ('U1', 'f8')})
        self.atomNames = atomData['names']
        self.atomMass  = atomData['mass']


    def parseCommonConfig(self, scf):
        """        
        Parses configuration for running molpro from the provided molpro config file 
        and sets up different methods and attributes relavant to the configuration 
        """
        # this just parses common configuration used in all the case. Systema and coordinate specific
        # configuration (like grid info) will be parsed in their respective inctances 
        molInfo =  dict(scf.items('molInfo'))
        
        self.scrdir = molInfo['scrdir']
        self.memory = molInfo['memory']
        try:
            self.proc = molInfo['processor']
        except KeyError:
            self.proc = '1'


        try:
            self.symmetry = scf.get('sysInfo', 'symmetry')
        except KeyError: # no symmetry keyword is provided, use `nosym`
            self.symmetry = 'nosym'
        self.nIREPs = self.IREP_LIST[self.symmetry]


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
            # if restricted keyword is provided then do `hf` else do an `uhf`
            if self.eInfo['restricted'].lower() == 'true':
                self.eInfo['hf'] = 'hf'
        except:
            pass
        if 'uhf_extra' in self.eInfo.keys():
            self.eInfo['hf'] += ';' + self.eInfo['uhf_extra']



    def parseResult(self, file):
        ''' Parses result from the ouptput .res files'''
        with open(file,"r") as f:
            dat = f.read().replace("D","E").strip().split("\n")[1:]
        dat = [[float(j) for j in i.strip().split()] for i in dat]
        return np.array(dat)




    def interp(self, file):
        ''' Fills the missing values in output file using a 1D interpolation '''
        # the file a 2D file
        data = np.loadtxt(file)
        grid1 = np.unique(data[:,0])   # numpy precision error
        grid2 = np.unique(data[:,1])



        nRows = grid2.shape[0]
        nCols = data.shape[1]
        res = np.array([], dtype=np.float64).reshape(0,nCols)

        for g1 in grid1:
            g1Block = data[data[:,0] == g1]

            if g1Block.shape[0]==1:   # g1=0 case, only one value is present. interpolation cant be applied here, copy the one single value everywhere
                res1 = np.full((nRows, nCols), g1Block)
                res1[:,1] == grid2 
            else:
                res1 = np.column_stack([np.full((nRows,), g1), grid2])
                for col in range(2, nCols):
                    gny  = self.interpolate(g1Block[:,1], g1Block[:,col], grid2)
                    res1 = np.column_stack([res1, gny])
            res = np.vstack([res, res1])
        return res



    def msg(self, m, cont=False):
        ''' Writes info in the log files'''
        if not cont : 
            m = '{:.<90}'.format(datetime.now().strftime("[%d-%m-%Y %I:%M:%S %p]     ") + m)
        else:
            m+='\n'
        self.logFile.write(m)
        self.logFile.flush()



    def writeFile(self, file, data):
        ''' Writes output data in plain txt'''
        file = open(file,'w')
        for tp in np.unique(data[:,0]):
            np.savetxt( file, data[data[:,0]==tp] ,delimiter="\t", fmt=str("%.8f"))
            file.write("\n")
        file.flush()



    def moveFiles(self, path):
        ''' Saves the geometry out and results file in a specific directory after running a job'''
        os.makedirs(path)
        files = glob('*.xyz') + glob('*.out') + glob('*.res') + glob('*.xml')
        for file in files:
            shutil.move(file, path)



    def cleanDirectory(self):
        # Clean and setup directory for new jobs
        #delete jobs folder if exists
        if os.path.isdir('IncompleteJobs'): shutil.rmtree('IncompleteJobs')
        if os.path.isdir('CompleteJobs')  : shutil.rmtree('CompleteJobs')

        # create jobs folder
        os.mkdir('IncompleteJobs')
        os.mkdir('CompleteJobs')


        # remove any existing files from other runs
        for pattern in ['*.res', '*.xyz', '*.xml*', '*.out*', '*.wfu', '*_irep_*']:
            for file in glob(pattern):
                try:
                    os.remove(file)
                except:
                    continue



    def runThisMolpro(self, grid1, gridn1, grid2, gridn2, filEe, fileN1, fileN2):
        '''Runs the molpro for each gridpoints '''

        # a general function to run the molpro, 
        # this function is agnostic about system and runs molpro for the grid passed through the grid1, grid2 arguments

        # subprocess calls blocks system I/O buffer, explicit I/O buffer counters are used to flush buffer
        self.logFile = open('adt_molpro.log', 'w')
        # open files to store result
        filee  = open(filEe, 'w', buffering=1)
        filen1 = open(fileN1,'w', buffering=1)
        filen2 = open(fileN2,'w', buffering=1)

        # Adding an info line at the top, so anyone can easily read the file about what colum is what
        # results for all the ireps are dumped into a single file and the header will keep the info about the irep
        # later it will be split into different ireps to do the adt
        eInfoLine  = nInfoLine = "#"+ gridn1.ljust(10) + ' \t ' + gridn2.ljust(10)
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                nInfoLine += '{a:>10}.{r}_{b}.{r}'.format(a=l,b=u,r=nIrep)
        nInfoLine += '\n'
        filen1.write(nInfoLine)
        filen2.write(nInfoLine)

        for nIrep, state in enumerate(self.nStateList, start=1):
            for s in range(1,state+1):
                eInfoLine += '\t{:>10}.{}'.format(s,nIrep)
        eInfoLine += '\n'
        filee.write(eInfoLine)



        # setup directory before run
        self.cleanDirectory()
        # run the initial point and save the wavefunction
        self.equiRun()

        done = False 

        for g2 in grid2:
            shutil.copy('molpro_init.wfu', 'molpro.wfu')      # copy the wfu from the initial/equilibrium job
            for g1 in grid1: 
                
                # if g1 = 0 (i.e. rho/theta/q) is provided then run g1=g2=0 only once
                # as for g1=0, all g2 will be same

                if (g1==0) & (g2!=0) & done:
                    continue # don't run the job for other rhos
                else:
                    done = True

                # this create geom function is modified inside respective inctances wrt system and coordinate type
                self.createGridGeom(g1, g2)
                # update the log keeping track of jobs
                self.msg('Running molpro job for {} = {}, {} = {}.......'.format(gridn1, g1, gridn2, g2))
                # path to save the job files
                path = 'CompleteJobs/{}_{}/{}_{}'.format(gridn1, g1, gridn2, g2)
                # copy the wfu file to restart the job
                shutil.copy('molpro.wfu', self.scrdir)
                # run the job if fails skip further of the loop
                exitcode = subprocess.call(
                    ['molpro', '-n', self.proc, "-d", self.scrdir, '-W .', 'grid.com']
                    )
                if exitcode==0:
                    self.msg( ' Job successful.', cont=True)
                else:
                    self.msg(' Job failed.', cont=True)
                    self.moveFiles(path.replace('C', "Inc"))
                    continue


                enrData = self.getEnergy()
                # this getTau function is created inside the respective class depending, system, coordinate and nact calculation type
                tau1, tau2 = self.getTau(g1, g2)
                # [None] increase one ndarray dimension so that savetxt can save it as row
                np.savetxt(filee, np.append([g1,g2],enrData)[None], fmt=str("%.8f"), delimiter='\t')
                np.savetxt(filen1, np.append([g1,g2],tau1)[None],  fmt=str("%.8f"), delimiter='\t')
                np.savetxt(filen2, np.append([g1,g2],tau2)[None],  fmt=str("%.8f"), delimiter='\t')
                self.moveFiles(path)
            filee.write('\n')
            filen1.write('\n')
            filen2.write('\n')
        # self.removeFiles(allOut=True)   # removes the wfu and .com files after complete run

        filee.close()
        filen1.close()
        filen2.close()

        # get the value to scale the energy data
        scalingVal = self.runScaleCalc()
        self.msg('All molpro jobs done.\n' ,cont=True)

        # indexes of columns to convert to radian
        # the first column will be turned to radian only when theta is there i.e scattering case for fixed rho
        # for all the other cases the first colum is an radial type coordinate, so not convert to radian
        # the attribute fixedRho exist only in scattering fixed rho case, so this check is enough
        # the second column is phi, so convert to radian anyway
        radCols = [0, 1] if  hasattr(self, 'fixedRho') else [1]



        # fill the missing values by 1D interpolation
        # '_mod' suffix means files with filled data
        inFiles = [filEe, fileN1, fileN2]

        for iFile in inFiles:
            dat = self.interp(iFile)
            # split after second column
            grid, data = np.split(dat, [2], axis=1)

            grid[:,radCols] = np.deg2rad(grid[:,radCols]) 

            # split the outputfiles according to their respective IRPEs
            # structure of enrgy file will be different that nacts, so handle that seperately
            # also scale the energy file data
            if iFile == filEe:    # energy file
                data -= scalingVal                   #     <<<<===== Scaling value depend on IREP
                ll = np.cumsum(self.nStateList)      # list containing the index where to split the columns wrt IREPs
            else:                 # nact files
                ll = np.cumsum(self.nTau)



            datas = np.split(data, ll, axis=1)
            for nIrep, dat in enumerate(datas, start=1):
                if not dat.size:        # balnk array, i.e this irep has no energy/nact calculation, skip
                    continue
                oFile = iFile.replace('.dat', '_irep_%s.dat'%nIrep)
                dat = np.column_stack([grid,dat])
                self.writeFile(oFile, dat)





    # decorator to scaling energy calculation, geometry is created in respective classes
    @classmethod
    def scaleWrapper(cls, func):
        def innerFunc(cls):
            if 'scale' not in cls.eInfo.keys():
                # Just returning a 0 for the non scaled calculation 
                # NOTE: in case of non sacled case only a single number is returned
                # while in sacled case, a array of number is returned 
                # but numpy handle both cases equaly during substraction.
                return 0
            func(cls)

            cls.msg( "Running molpro job for scaling geometry.......")
            exitcode = subprocess.call(
                ['molpro', '-n', cls.proc, "-d", cls.scrdir, '-W .', 'scale.com']
                )
            if exitcode==0: 
                cls.msg( ' Job successful.', cont=True)
                scale = cls.getEnergy()
                cls.moveFiles('CompleteJobs/Scaling_point')
                # now the scale has to be consistent with the IREPs, 
                # so make the scale a 1D array with same IREPs having the same value, i.e. lowest state value
                ll = np.split(scale, np.cumsum(cls.nStateList)[:-1])  #split energies in list of different IREPs
                for i in ll:
                    if not i.size: continue
                    i[:] = i[0]                                        # subtract IREP from ground state
                scale = np.append(*ll)                                 # stitch it up
                return scale
            else:
                cls.msg( ' Job failed', cont=True)
                cls.moveFiles('IncompleteJobs/Scaling_point')
                return 'not done'
        return innerFunc



    # decorator to initial point calculation, geometry is created in respective classes
    @classmethod
    def initWrapper(cls, func):
        def innerFunc(cls):
            func(cls)
            exitcode = subprocess.call(
                ['molpro', '-n', cls.proc, "-d", cls.scrdir, '-W .', 'init.com']
                )
            if exitcode==0: 
                cls.msg( ' Job successful.', cont=True)
                cls.moveFiles('CompleteJobs/Equilibrium_point')
            else:
                cls.msg( ' Job failed', cont=True)
                sys.exit('Molpro failed in equilibrium step')
        return innerFunc




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
        assert len(self.rhoList)==3, "A grid of rho is required for scattering case"
        self.rhoGrid = self.makeGrid(self.rhoList)

        self.phiList = [float(i) for i in self.gInfo['phi'].split(',')]
        assert len(self.phiList)==3, "A grid of phi is required for scattering case"
        self.phiGrid = self.makeGrid(self.phiList)

        # setup proper function for nact type
        if self.nInfo['method']=='cpmcscf':
            self.createGridGeom = self.createOneGeom
            self.getTau         = self.getTauAna


        elif self.nInfo['method'] == 'ddr':
            try:
                self.d1             = float(self.gInfo['drho'])
                self.d2 = float(self.gInfo['dphi'])
            except KeyError as key:
                raise Exception("%s keyword as geometry increment is required for ddr NACT calculation"%str(key))

            self.createGridGeom = self.createAllGeom
            self.getTau         = self.getTauDdr

        self.createTemplate()


    def parseSData(self, geomFile, freqFile, wilsonFile):
        '''Parses equilibrium geometry, frequencies and the wilson matrix data for a sceptroscopic system'''
        self.equiGeom = np.loadtxt(geomFile)
        freq = np.loadtxt(freqFile)
        wilson = np.loadtxt(wilsonFile)

        wilson = wilson.reshape(self.equiGeom.shape[0], 3, freq.shape[0])
        freqInv = np.sqrt(self.hbar/(freq*self.cInvToTauInv))
        massInv = np.sqrt(1 / self.atomMass)
        # an array containing information about wilson, frequency and mass, used for geometry creation and parsing analytic nact
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
        self.createOneGeom(rho,  phi,  'geom.xyz')
        self.createOneGeom(rho+self.d1, phi, 'geom2.xyz')
        self.createOneGeom(rho-self.d1, phi, 'geom3.xyz')
        self.createOneGeom(rho, phi+self.d2,'geom4.xyz')
        self.createOneGeom(rho, phi-self.d2,'geom5.xyz')



    def getTauAna(self,rho,phi):
        # So this will ultimately return two list one for taurho and another for tauphi
        # each list will be ordered as 12,13,23...(IREP1),12,13,23...(IREP2) etc. 
        # the same order in the list self.nactPairList will be maintaned everywhere
        rotoMat = np.array([[self.cos(phi), self.sin(phi)], [-rho*self.sin(phi), rho*self.cos(phi)]])
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ananac_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                grads = self.angtobohr*self.parseResult(file)
                tau = np.einsum('ijk,ij,lk->l',self.wilFM[...,self.vModes], grads, rotoMat)
                tauList.append(tau)
        return np.abs(tauList).T



    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ddrnact_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                tau = self.parseResult(file)
                tauList.append(tau)
        return np.abs(tauList).T



    @Base.initWrapper
    def equiRun(self):
        '''Runs molpro for the equilibrium geometry'''
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program for equilibrium geometry.\n"
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *self.equiGeom[i])

        with open('geom.xyz',"w") as f: f.write(tmp)

        self.msg( "Running molpro job for equilibrium point.......")



    @Base.scaleWrapper
    def runScaleCalc(self):
        scaleList = [float(i) for i in self.eInfo['scale'].split(',')]
        assert len(scaleList)==2, "Provide coordinate for scale energy calculation"
        rho, phi = scaleList
        self.createOneGeom(rho,phi)



    def runMolpro(self):
        file = ['energy.dat','tau_rho.dat','tau_phi.dat']
        self.runThisMolpro(
            self.rhoGrid, 
            'Rho', 
            self.phiGrid, 
            'Phi', *file )

        # len(nTau) = len(nstate) = number of IREPs
        fileList=[]
        for nIrep, (enr, tau) in enumerate(zip(self.nTau, self.nStateList), start=1):
            tmp = []
            if enr:   # If enr 0, means no energy is calculated for this IREP
                tmp.append( "energy_irep_{}.dat".format(nIrep) )
            if tau:   # If enr 0, means no NACT is calculated for this IREP
                tmp.extend(
                   ["tau_rho_irep_{}.dat".format(nIrep) ,
                    "tau_phi_irep_{}.dat".format(nIrep) ]
                )
            # when state is 0, this actually append an empty list
            fileList.append(tmp)
        return ("2D", fileList)



class Scattering(Base):
    ''' Inherited from the Base class, this class containg necessary methods
     for running molpro for a Scattering system in hyperspherical coordinate'''
    def __init__(self, conFig, atomFile):
        self.parseData(atomFile)
        self.parseCommonConfig(conFig)
        self.parseThisConfig(conFig)


    def parseThisConfig(self, conFig):
        # implement rho, theta phi all as variable.
        rhoList = [float(i) for i in self.gInfo['rho'].split(',')]
        thetaList = [float(i) for i in self.gInfo['theta'].split(',')]
        r= len(rhoList)
        t= len(thetaList)

        if (r==1) & (t==3): # ab initio will be done for a fixed rho
            self.fixedRho = True
            self.rho = rhoList[0]
            self.thetaGrid = self.makeGrid(thetaList)
        elif (r==3) & (t==1): # ab initio will be done for a fixed theta
            self.fixedRho = False
            self.theta = thetaList[0]
            self.rhoGrid =  self.makeGrid(rhoList)
        else:
            raise Exception("Provide a grid for either Rho or Theta and give a fixed value of the other for the scattering case")


        phiList = [float(i) for i in self.gInfo['phi'].split(',')]
        assert len(phiList)==3, "A grid of phi is required for the scattering case"
        self.phiGrid = self.makeGrid(phiList)

        if self.nInfo['method']=='cpmcscf':
            self.getTau = self.getTauAna
            self.createGridGeom = self.createOneGeom

        elif self.nInfo['method']=='ddr':
            try:
                self.d1 = float(self.gInfo['dtheta']) if self.fixedRho else float(self.gInfo['drho'])
                self.d2 = float(self.gInfo['dphi'])
            except KeyError as key:
                raise Exception("%s keyword as geometry increment is required for ddr NACT calculation"%str(key))
            self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr
    
        self.createTemplate()


    def AreaTriangle(self,a,b,c):
        """ area of a tringle with sides a,b,c """
        ps = (a+b+c)/2.0
        ar = ps*(ps-a)*(ps-b)*(ps-c)
        # negative area due to round off errors set to zero
        if ar < 0.0:
            ar = 0.0
        ar = np.sqrt(ar)
        return ar

    def toJacobi(self,rho, theta,phi):
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
        R1 = (1.0/np.sqrt(2.0))*rho*d3*np.sqrt(1.0+ self.sin(theta)*self.cos(phi+eps3)) 
        R2 = (1.0/np.sqrt(2.0))*rho*d1*np.sqrt(1.0+ self.sin(theta)*self.cos(phi))      
        R3 = (1.0/np.sqrt(2.0))*rho*d2*np.sqrt(1.0+ self.sin(theta)*self.cos(phi-eps2)) 

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

    def hyperToCart(self, rho, theta, phi):
        # create an cartesian corodinate from hyper spherical coordiante
        rs, rc, gamma = self.toJacobi(rho, theta, phi)
        p1 = [0, rc*np.sin(gamma), rc*np.cos(gamma)]
        p2 = [0,0.0, -rs/2.0]
        p3 = [0,0.0, rs/2.0 ]
        return np.array([p1, p2, p3])*self.bohrtoang # return in angstrom


    # The argument `rhoTheta` is actually for passing both rho or theta values to this function
    # depending on the situation and making it compatible for both varying rho and theta
    # Why this complicacy you ask? to make this consistent with the geometry creation function
    # from other objects as they all take two variable and making it a bit more concise

    def createOneGeom(self, rhoTheta, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given theta , phi'''
        if self.fixedRho:
            rho = self.rho
            theta = rhoTheta
        else:
            rho = rhoTheta
            theta = self.theta

        curGeom  = self.hyperToCart(rho, theta, phi)
        msg = 'for Rho = {}, Theta={}, Phi = {}'.format(rho, theta, phi)
        nAtoms = len(self.atomNames)
        tmp = " {}\n".format(nAtoms)
        tmp+= "Geometry file created from ADT program. %s \n"%msg
        for i in range(nAtoms):
            tmp += "{},{},{},{}\n".format(self.atomNames[i], *curGeom[i])
        with open(outFile,"w") as f:
            f.write(tmp)


    # following function is used when fixed theta and rho is varying 
    def createAllGeom(self, rhoTheta, phi):
        ''' Creates 5 different geometry files for using in molpro ddr calculation '''
        self.createOneGeom(rhoTheta,    phi,    'geom.xyz')
        self.createOneGeom(rhoTheta+self.d1, phi,   'geom2.xyz')
        self.createOneGeom(rhoTheta-self.d2, phi,   'geom3.xyz')
        self.createOneGeom(rhoTheta,    phi+self.d2,'geom4.xyz')
        self.createOneGeom(rhoTheta,    phi-self.d2,'geom5.xyz')



    def getTauAna(self, rhoTheta, phi):
        '''Used in Analytical NACT calculation'''

        # the analytical nact is calculated by numercally differentiating the grid with an increment of 1/100 of grid spacing
        #  as cartesian to hepersphercial conversion is ? 
        if self.fixedRho:
            rho = self.rho
            theta = rhoTheta
            dTheta = (self.thetaGrid[1]-self.thetaGrid[0])/100.0
            gradThetaPlus  = self.hyperToCart(rho, theta+dTheta, phi)
            gradThetaMinus = self.hyperToCart(rho, theta-dTheta, phi)
            gradRhoTheta      = (gradThetaPlus - gradThetaMinus)/2*dTheta
        else:
            rho = rhoTheta
            theta = self.theta
            dRho = (self.rhoGrid[1]-self.rhoGrid[0])/100.0
            gradRhoPlus  = self.hyperToCart(rho+dRho,theta, phi)
            gradRhoMinus = self.hyperToCart(rho-dRho,theta, phi)
            gradRhoTheta      = (gradRhoPlus - gradRhoMinus)/2*dRho

        dPhi = (self.phiGrid[1]-self.phiGrid[0])/100.0
        gradPhiPlus  = self.hyperToCart(rho, theta, phi+dPhi)
        gradPhiMinus = self.hyperToCart(rho, theta, phi-dPhi)
        gradPhi      = (gradPhiPlus - gradPhiMinus)/2*dPhi 



        # rotoMat = np.array([[self.cos(phi), self.sin(phi)], [-rho*self.sin(phi), rho*self.cos(phi)]])
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ananac_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                grads = self.angtobohr*self.parseResult(file)
                tauRhoTheta = np.einsum('ij,ij', grads, gradRhoTheta)
                tauPhi   = np.einsum('ij,ij', grads, gradPhi)
                tauList.append([tauRhoTheta, tauPhi])
        return np.abs(tauList).T



    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ddrnact_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                tau = self.parseResult(file)
                tauList.append(tau)
        return np.abs(tauList).T


    @Base.initWrapper
    def equiRun(self):
        '''Runs molpro for a initial geometry, i.e. theta=phi=0'''
        if self.fixedRho:
            self.createOneGeom(self.thetaGrid[0], self.phiGrid[0])
        else:
            self.createOneGeom(self.rhoGrid[0], self.phiGrid[0])
        self.msg( "Running molpro job for initial point....." )




    @Base.scaleWrapper
    def runScaleCalc(self):
        scaleList = [float(i) for i in self.eInfo['scale'].split(',')]
        assert len(scaleList)==3, "Provide coordinate for scale energy calculation"
        rho,theta, phi = scaleList
        if self.fixedRho:
            self.rho =rho
            self.createOneGeom(theta,phi)
        else:
            self.theta = theta
            self.createOneGeom(rho,phi)
        



    def runMolpro(self):
        if self.fixedRho:
            file = ['energy.dat','tau_theta.dat','tau_phi.dat']
            self.runThisMolpro(
                self.thetaGrid, 
                'Theta', 
                self.phiGrid, 
                'Phi', *file )
            tXt="theta"
        else :
            file = ['energy.dat','tau_rho.dat','tau_phi.dat']
            self.runThisMolpro(
                self.rhoGrid, 
                'Rho', 
                self.phiGrid, 
                'Phi', *file )
            tXt = 'phi'

        # now len(nTau) = len(nstate) = number of IREPs
        fileList=[]
        for nIrep, (enr, tau) in enumerate(zip(self.nTau, self.nStateList), start=1):
            tmp = []
            if enr:   # If enr 0, means no energy is calculated for this IREP
                tmp.append( "energy_irep_{}.dat".format(nIrep) )
            if tau:   # If enr 0, means no NACT is calculated for this IREP
                tmp +=  ["tau_{}_irep_{}.dat".format(tXt, nIrep), "tau_phi_irep_{}.dat".format(nIrep) ]
            fileList.append(tmp)
        return ("2D", fileList)




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
        else :            # 2D Jacobi will be done
            self.Jacobi1D = False
            ql = [float(i) for i in ql]
            self.qGrid = self.makeGrid(ql)



        if self.nInfo['method']=='cpmcscf':
            self.getTau = self.getTauAna
            if self.Jacobi1D:
                self.createGridGeom = self.createOneGeom1D
            else:
                self.createGridGeom = self.createOneGeom

        elif self.nInfo['method']=='ddr':
            if self.Jacobi1D:
                try:
                    self.d1 = float(self.gInfo['dphi'])
                except KeyError as key:
                    raise Exception("%s keyword as geometry increment is required for ddr NACT calculation"%str(key))
                self.createTemplate = self.createDdrTemplate1D
                self.createGridGeom = self.createAllGeom1D
            else:
                try:
                    self.d1 = float(self.gInfo['dq'])
                    self.d2 = float(self.gInfo['dphi'])
                except KeyError as key:
                    raise Exception("%s keyword as geometry increment is required for ddr NACT calculation"%str(key))
                self.createGridGeom = self.createAllGeom
            self.getTau = self.getTauDdr

        self.createTemplate()



    def createDdrTemplate1D(self):
        # this ddr template is  done by a creating a 2D ddr template and modifying it
        self.d2 = self.d1
        # self.createTemplate called above is actually refers to this function for 1D ddr case 
        # The following calls the createTemplate from the parent class, which actually creates ddr template for 2D case
        super(Jacobi, self).createTemplate()   #! <<<<==== CAUTION: ordering of the inheritence methods

        with open("grid.com","r") as f:
            template = f.read()
        # remove 7 lines starting from !N1D and also remove the word ",nacmp"
        # regular expression on the rescue
        odTemplate = re.sub("!N1D:((.*\n){7})|,nacmp", '', template)
        odTemplate = re.sub("!N1DB.*!N1DB", '', odTemplate, flags=re.DOTALL|re.MULTILINE)

        with open("grid.com","w") as f:
            f.write(odTemplate)



    def getTauAna(self, *args):
        '''Used in Analytical NACT calculation'''
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ananac_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                grads = self.parseResult(file)
                if self.Jacobi1D:
                    phi = args[0] # 1D case only phi is passed
                    tau = self.q*(-grads[2,0]*self.sin(phi) + grads[2,1]*self.cos(phi))
                    tauList.append(tau)
                else :
                    q, phi = args # 2D case q,phi is passed
                    tauphi = q*(-grads[2,0]*self.sin(phi) + grads[2,1]*self.cos(phi))
                    tauq   = grads[2,0]*self.cos(phi) + grads[2,1]*self.sin(phi)
                    tauList.append([tauphi, tauq])
        return np.abs(tauList).T   #! CAUTION about transposing in 1D case


    def getTauDdr(self, *args):
        '''Used in DDR NACT calculation'''
        # args are not really used here, but to keep consistency
        tauList = []
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                file = "ddrnact_{a}.{r}_{b}.{r}.res".format(a=l,b=u,r=nIrep)
                tau = self.parseResult(file)
                tauList.append(tau)
        return np.abs(tauList).T


    def createOneGeom(self, q, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given r, R, gamma, q, phi'''

        # three particle is on yz plane with the diatom being on the z axis
        curGeom  = self.bohrtoang*np.array([
            [0.0, 0.0, -self.smallR/2.0],
            [0.0, 0.0, self.smallR/2.0],
            [0.0, self.capR*self.sin(self.gamma)+q*self.sin(phi), self.capR*self.cos(self.gamma)+q*self.cos(phi)]])


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
        self.createOneGeom(q, phi, 'geom.xyz')
        self.createOneGeom(q+self.d1, phi, 'geom2.xyz')
        self.createOneGeom(q-self.d1, phi, 'geom3.xyz')
        self.createOneGeom(q, phi+self.d2,'geom4.xyz')
        self.createOneGeom(q, phi-self.d2,'geom5.xyz')



    def createOneGeom1D(self, phi, outFile='geom.xyz'):
        ''' Creates geometry file, to be used in molpro for the given r, R, gamma, q, phi'''
        # same oneGeomTemplate but with a fixed q
        self.createOneGeom(self.q, phi, outFile)



    def createAllGeom1D(self, phi):
        ''' Creates 3 different geometry files for using in molpro ddr calculation '''

        self.createOneGeom1D(phi,    'geom.xyz')
        self.createOneGeom1D(phi+self.d1,'geom2.xyz')
        self.createOneGeom1D(phi-self.d1,'geom3.xyz')



    @Base.scaleWrapper
    def runScaleCalc(self):
        scaleList = [float(i) for i in self.eInfo['scale'].split(',')]
        assert len(scaleList)==3, "Provide coordinate for scale energy calculation"
        self.smallR, self.capR, self.gamma = scaleList
        # scaling is run at the end so modifying r,R,gamma globally and passing rho=phi=0
        if self.Jacobi1D:
            self.createOneGeom(0)
        else:
            self.createOneGeom(0,0)


    # only called for 2D Jacobi study
    @Base.initWrapper
    def equiRun(self):
        '''Runs molpro for a initial geometry, i.e. theta=phi=0'''
        self.createOneGeom(self.qGrid[0],self.phiGrid[0])
        self.msg( "Running molpro job for initial point....." )


    def runMolpro(self):
        file = ['energy.dat','tau_q.dat','tau_phi.dat']
        self.runThisMolpro(
            self.qGrid, 
            'q', 
            self.phiGrid, 
            'Phi', *file )
        # now len(nTau) = len(nstate) = number of IREPs
        fileList=[]
        for nIrep, (enr, tau) in enumerate(zip(self.nTau, self.nStateList), start=1):
            tmp = []
            if enr:   # If enr 0, means no energy is calculated for this IREP
                tmp.append( "energy_irep_{}.dat".format(nIrep) )
            if tau:   # If enr 0, means no NACT is calculated for this IREP
                tmp.extend(
                    ["tau_q_irep_{}.dat".format(nIrep), "tau_phi_irep_{}.dat".format(nIrep) ]
                )
            fileList.append(tmp)
        return ("2D", fileList)


    def runMolpro1D(self):
        '''Runs the molpro for each gridpoints '''
        # this here is only for 1d jacobi
        # subprocess calls blocks system I/O buffer, so the I/Os (sometimes) have to be flushed out explicitely 


        self.logFile = open('adt_molpro.log', 'w')

        # setup directory before run
        self.cleanDirectory()


        # open files to store result
        filee  = open('energy.dat', 'w', buffering=1)
        filen  = open('tau_phi.dat','w', buffering=1)


        # Adding an info line at the top, so anyone can easily read the file
        eInfoLine  = nInfoLine = "#          Phi"
        for nIrep, nactPairs in enumerate(self.nactPairsList, start=1):
            for l,u in nactPairs:
                nInfoLine += '{a:>10}.{r}_{b}.{r}'.format(a=l,b=u,r=nIrep)
        nInfoLine += '\n'
        filen.write(nInfoLine)

        for nIrep, state in enumerate(self.nStateList, start=1):
            for s in range(1,state+1):
                eInfoLine += '\t{:>10}.{}'.format(s,nIrep)
        eInfoLine += '\n'
        filee.write(eInfoLine)
        # firsDone True means one job is done so, copy wfu file to restart the new job
        firstDone = False



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
                self.moveFiles(path.replace('C', "Inc"))
                continue 
            firstJobDone = True
            enrData = self.getEnergy()
            tau = self.getTau(phi)
            self.moveFiles(path)
            np.savetxt(filee,  np.append([phi],enrData)[None], fmt=str("%.8f"), delimiter='\t')
            np.savetxt(filen, np.append([phi],tau)[None],  fmt=str("%.8f"), delimiter='\t')
        
        scalingVal = self.runScaleCalc()

        self.msg('All molpro jobs done.', cont=True)

        inFiles = ['energy.dat', 'tau_phi.dat']

        for iFile in inFiles:
            dat = np.loadtxt(iFile)
            grid = self.phiGrid
            data = np.column_stack([self.interpolate(dat[:,0], dat[:,col], grid) for col in range(1, dat.shape[1])])

            grid = np.deg2rad(grid)        # second column is radian, convert it to phi

            # split the outputfiles according to their respective IRPEs
            # structure of enrgy file will be different that nacts, so handle that seperately
            # also scale the energy file data
            if iFile == 'energy.dat':    # energy file
                data -= scalingVal                   #     <<<<===== Scaling value depend on IREP
                ll = np.cumsum(self.nStateList)      # list containing the index where to split the columns wrt IREPs
            else:                 # nact files
                ll = np.cumsum(self.nTau)


            datas = np.split(data, ll[:-1], axis=1)
            for nIrep, dat in enumerate(datas, start=1):
                oFile = iFile.replace('.dat', '_irep_%s.dat'%nIrep)
                dat = np.column_stack([grid,dat])
                np.savetxt(oFile, dat, fmt=str("%.8f"), delimiter='\t')

        fileList=[]
        for nIrep, (enr, tau) in enumerate(zip(self.nTau, self.nStateList), start=1):
            tmp = []
            if enr:   # If enr 0, means no energy is calculated for this IREP
                tmp.append( "energy_irep_{}.dat".format(nIrep) )
            if tau:   # If enr 0, means no NACT is calculated for this IREP
                tmp.append("tau_phi_irep_{}.dat".format(nIrep) )
            fileList.append(tmp)
        return ("1D", fileList)


