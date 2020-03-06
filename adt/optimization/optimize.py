from __future__ import absolute_import, unicode_literals, division, print_function
import re
import sys
import textwrap
import subprocess
import numpy as np

if sys.version_info.major>2:
    from configparser import ConfigParser as ConfigParser
else :
    from ConfigParser import SafeConfigParser as ConfigParser

class GaussianOptg():

    def __init__(self, config):
        scf = ConfigParser({'path':'g16','symmetry':'nosym','processor':'1'})
        scf.read(config)

        self.optInfo  = dict(scf.items('optInfo'))
        for key in ['memory', 'basis', 'method', 'spin', 'charge', 'symmetry']:
            if key not in self.optInfo:
                raise Exception('Option "%s" not found in config'%key)
        try:
            self.geomFile = scf.get('gInfo','file')
        except:
            raise KeyError('Initial geometry file not found in config')
        self.CreateTemplate()

    def CreateTemplate(self):
        gaussianTemplate = textwrap.dedent('''            %nprocshared={proc}
            %chk=optg.chk
            %Mem={memory}
            #opt=(cartesian,maxcycles=50) freq {method}/{basis} scf=qc {sym}

            title name

            {charge} {spin}
            '''.format(proc   = self.optInfo['processor'],
                       memory = self.optInfo['memory'],
                       method = self.optInfo['method'],
                       basis  = self.optInfo['basis'],
                       sym    = self.optInfo['symmetry'],
                       spin   = self.optInfo['spin'],
                       charge = self.optInfo['charge']))

        with open(self.geomFile) as f: geomDat = f.read()
        gaussianTemplate += geomDat + '\n'
        with open('optg.inp', 'w') as f: f.write(gaussianTemplate)

    def runOpt(self):
        try:
            rungaussian = subprocess.call([self.optInfo['path'],'optg.inp'])
        except:
            raise Exception('Can not run Gaussian.')
        if rungaussian:
            raise Exception(' Optimization failed.')

    def getResults(self):
        with open('optg.log') as f: txt = f.read()
        atomNo = int(re.findall(r'NAtoms=\s+(\d+)',txt)[0])
        optGeom = re.findall(r' Optimization completed.(?:(?:.*\n){{9}})((?:.*\n){{{}}})'.format(atomNo), txt)
        optGeom = [i.split() for i in filter(None,optGeom[0].split('\n'))]
        optGeom = np.array(optGeom,dtype=np.float64)[:,3:]
        np.savetxt('equigeom.dat', optGeom, fmt=str('%.10f'), delimiter='\t')

        freDat= re.findall(r' Frequencies --(.*)(?:(?:.*\n){{5}})((?:.*\n){{{}}})'.format(atomNo),txt)
        freq = np.array( [j for i,_ in freDat for j in i.split() ], dtype=np.float)
        # freq = np.array( [i.split() for i,_ in freDat], dtype=np.float64).reshape(-1)
        np.savetxt('frequency.dat', freq[None,:], fmt=str('%.10f'), delimiter='\t')

        # # remove columns regarding the atom number and atom type (one row for one atom)
        tmp= [[i.split()[2:] for i in filter(None,j.split('\n')) ] for _,j in freDat ]
        wilson = np.asarray(np.column_stack(tmp), dtype=np.float)
        # structure the array, as required by the ADT code
        res = np.vstack(np.transpose(wilson.reshape(atomNo,-1,3),(0,2,1)))
        np.savetxt('wilson.dat', res, fmt=str('%.10f'), delimiter='\t')


class GamessOptg():

    def __init__(self, config):
        scf = ConfigParser({'path':'rungms','symmetry':'c1','processor':'1'})
        scf.read(config)

        self.optInfo  = dict(scf.items('optInfo'))
        for key in ['memory', 'memddi', 'basis', 'method', 'spin', 'charge', 'symmetry']:
            if key not in self.optInfo:
                raise Exception('Option "%s" not found in config'%key)
        try:
            self.geomFile = scf.get('gInfo','file')
        except:
            raise KeyError('Initial geometry file not found in config')
        self.CreateTemplate()

    def CreateTemplate(self):
        indent = lambda txt:'\n'.join([' '+i.strip() for i in filter(None,txt.split('\n'))])
        with open(self.geomFile) as f: geomDat = f.read()
        self.nAtoms = len(filter(None,geomDat.split('\n')))

        gamessTemplate = textwrap.dedent('''
                                            $CONTRL SCFTYP={scfmeth} {lvl} RUNTYP=OPTIMIZE ICHARG={charge}
                                            COORD=UNIQUE MULT={spin} MAXIT=200 ISPHER=1 $END
                                            $SYSTEM MWORDS={memory} MEMDDI ={memddi} $END
                                            {pre}
                                            $STATPT NSTEP=100 HSSEND=.T. $END
                                            {post}
                                            $BASIS  GBASIS={basis} $END
                                            $GUESS  GUESS=HUCKEL $END
                                            $DATA
                                            optg and freq
                                            C1
                                            {geom}
                                            $END
                                            '''.format( scfmeth = 'RHF' if self.optInfo['spin'] =='1' else
                                                                  'UHF' if self.optInfo['method'] =='b3lyp' else 'ROHF',
                                                        lvl     = {'mp2':'MPLEVL=2','ump2':'MPLEVL=2',
                                                                  'ccsd':'CCTYP=ccsd','uccsd':'CCTYP=ccsd',
                                                                  'b3lyp':'DFTTYP=b3lyp'}.get(self.optInfo['method']),
                                                        charge  = self.optInfo['charge'],
                                                        spin    = self.optInfo['spin'],
                                                        memory  = self.optInfo['memory'],
                                                        memddi  = self.optInfo['memddi'],
                                                        pre     = '$SCF DIRSCF=.TRUE. $END\n$CPHF CPHF=AO $END'
                                                                  if self.optInfo['method'] == 'b3lyp' else '',
                                                        post    = '$CCINP MAXCC =100 $END\n$FORCE METHOD=FULLNUM $END'
                                                                  if self.optInfo['method'] in ['ccsd','uccsd'] else '',
                                                        basis   = self.optInfo['basis'],
                                                        geom    = geomDat.strip()))

        with open('optg.inp', 'w') as f: f.write(indent(gamessTemplate))

    def runOpt(self):
        try:
            rungamess = subprocess.call('{} optg.inp >& optg.log'.format(self.optInfo['path']), shell = True)
        except:
            raise Exception('Can not run gamess')
        if rungamess:
            raise Exception(' Optimization failed.')

    def getResults(self):
        with open('./optg.log') as f : txt = f.read()
        ind= int((re.findall(r' MODE\(?S\)?\s+\d+ TO\s+(\d+) ARE TAKEN AS ROTATIONS AND TRANSLATIONS.',txt))[0])
        optGeom = re.findall(r' EQUILIBRIUM GEOMETRY LOCATED.(?:(?:.*\n){{4}})((?:.*\n){{{}}})'.format(self.nAtoms), txt)
        optGeom = [i.split()[2:] for i in filter(None,optGeom[0].split('\n'))]
        optGeom = np.array(optGeom,dtype=np.float64)
        np.savetxt('equigeom.dat', optGeom, fmt=str('%.10f'), delimiter='\t')
        freq = [ float(j)  for i in re.findall('FREQUENCY:(.*)',txt) for j in i.replace('I','').split() ][ind:]
        freq = np.asarray(freq,dtype=np.float)
        np.savetxt('frequency.dat', freq[None,:], fmt=str('%.10f'), delimiter='\t')
        pat1 = r'FREQUENCY:\s+(.*)(?:(?:.*\n){{5}})((?:.*\n){{{}}})'.format(3*self.nAtoms) # 3*nAtoms
        pat2 = r'(?:\s+\d+\s+\w+)?\s+[X|Y|Z](.*)'
        freDat = re.findall(pat1,txt)
        arr =np.column_stack([ np.vstack([j.split() for j in re.findall(pat2,i)]) for _,i in freDat])
        wilson = np.array(arr, dtype=np.float)[:,ind:]
        np.savetxt('wilson.dat', wilson, fmt=str('%.10f'), delimiter='\t')

class MolproOptg(object):

    def __init__(self, config):
        scf = ConfigParser({'path':'molpro','symmetry':'nosym','processor':'1','charge':'0'})
        scf.read(config)

        self.optInfo  = dict(scf.items('optInfo'))
        for key in ['memory', 'basis', 'method', 'symmetry']:
            if key not in self.optInfo:
                raise Exception('Option "%s" not found in config'%key)
        try:
            self.geomFile = scf.get('gInfo','file')
        except:
            raise Exception('Initial geometry file not found in config')
        self.CreateTemplate()

    def CreateTemplate(self):
        method = self.optInfo['method']
        meth = {'mp2':'hf','ump2':'uhf','uccsd':'rhf','ccsd':'hf'}.get(method,'')
        method = ('ks,' if method=='b3lyp' else '') + method

        molproTemplate = textwrap.dedent('''
            memory,{memory}

            set,charge={charge}
            basis={basis}

            symmetry,{sym}

            geometry=geom.xyz
            {exmeth}
            {{{method}}}
            optg;
            {{frequencies;print,hessian}}

            ---
            '''.format(exmeth = meth,
                       charge = self.optInfo['charge'],
                       memory = self.optInfo['memory'],
                       basis  = self.optInfo['basis'],
                       sym    = self.optInfo['symmetry'],
                       method = method))

        with open(self.geomFile) as f: geomDat = f.read()
        self.nAtoms = len(filter(None,geomDat.split('\n')))

        with open('geom.xyz','w') as f:
            f.write(
                "{}\ngeometry file created for optimization\n{}".format(
                        self.nAtoms,
                        re.sub('[ \t]+',',',geomDat)
                        )
                    )

        with open('optg_mol.inp', 'w') as f: f.write(molproTemplate)

    def runOpt(self):
        try:
            runMolpro = subprocess.call([self.optInfo['path'], '-n',self.optInfo['processor'], 'optg_mol.inp','--no-xml-output'])
        except:
            raise Exception('Can not run Molpro')
        if runMolpro:
            raise Exception(' Optimization failed.')

    def getResults(self):
        with open('./optg_mol.out') as f : txt = f.read()

        nAtoms = self.nAtoms # int(re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){4})\s+(\d+)', txt)[0])
        optGeom = re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){{6}})((?:.*\n){{{}}})'.format(nAtoms), txt)[0]

        optGeom = [i.split()[1:] for i in filter(None, optGeom.split('\n'))]
        optGeom = np.array(optGeom, dtype=np.float)
        np.savetxt('equigeom.dat', optGeom, fmt=str('%.10f'), delimiter='\t')

        freDat= re.findall(r' Wavenumbers \[cm-1\](.*)(?:(?:.*\n){{3}})((?:.*\n){{{}}})'.format(3*nAtoms),txt)

        freqs = np.array([j for i,_ in freDat for j in i.split() ], dtype=np.float)
        freqs = freqs[freqs!=0] # take the vibrational freqs only i.e. non zero
        np.savetxt('frequency.dat', freqs[None,:], fmt=str('%.10f'), delimiter='\t')

        tmp= [[i.split()[1:] for i in filter(None,j.split('\n')) ] for _,j in freDat ]
        wilson = np.asarray(np.column_stack(tmp), dtype=np.float)[:,:freqs.shape[0]]
        np.savetxt('wilson.dat', wilson, fmt=str('%.10f'), delimiter='\t')





if __name__ == "__main__":
    p = GamessOptg('gms.config')
