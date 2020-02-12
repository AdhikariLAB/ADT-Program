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
        scf = ConfigParser({'symmetry':'nosym','processor':'1'})
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
        GaussianTemplate = textwrap.dedent('''            %nprocshared={proc}
            %chk=optg.chk
            %Mem={memory}
            #opt=cartesian freq {method}/{basis} scf=qc {sym}

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
        GaussianTemplate += geomDat + '\n'
        with open('optg.inp', 'w') as f: f.write(GaussianTemplate)
    
    def runOpt(self):
        try:
            rungaussian = subprocess.call(['g16','optg.inp'])
        except:
            raise Exception('Can not run Gaussian 16')
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




class MolproOptg(object):

    def __init__(self, config):
        scf = ConfigParser({'symmetry':'nosym','processor':'1'})
        scf.read(config)

        self.optInfo  = dict(scf.items('optInfo'))
        for key in ['memory', 'basis', 'method', 'symmetry']:
            if key not in self.optInfo:
                raise Exception('Option "%s" not found in config'%key)
        try:
            self.geomFile = scf.get('gInfo','file')
        except:
            raise KeyError('Initial geometry file not found in config')
        self.CreateTemplate()

    def CreateTemplate(self):
        GaussianTemplate = textwrap.dedent('''            
            memory,{memory}

            basis={basis}

            symmetry,{sym}

            geometry=geom.xyz
            
            {{{method}}}
            optg;
            {{frequencies;print,hessian}}

            ---
            '''.format(memory = self.optInfo['memory'],
                       basis  = self.optInfo['basis'],
                       sym    = self.optInfo['symmetry'],
                       method = self.optInfo['method']))

        with open(self.geomFile) as f: geomDat = f.read()
        self.nAtoms = len(filter(None,geomDat.split('\n')))

        with open('geom.xyz','w') as f:
            f.write(
                "{}\ngeometry file created for optimization\n{}".format(
                        self.nAtoms, 
                        re.sub('[ \t]+',',',geomDat)
                        )
                    )

        with open('optg_mol.inp', 'w') as f: f.write(GaussianTemplate)

    def runOpt(self):
        try:
            runMolpro = subprocess.call(['molpro', '-n',self.optInfo['processor'], 'optg_mol.inp','--no-xml-output'])
        except:
            raise Exception('Can not run Molpro')
        if runMolpro:
            raise Exception(' Optimization failed.')

    def getResults(self):
        with open('./optg_mol.out') as f : txt = f.read()

        nAtoms = self.nAtoms # int(re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){4})\s+(\d+)', txt)[0])
        optGeom,= re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){{6}})((?:.*\n){{{}}})'.format(nAtoms), txt)

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
    p = GaussianOptg('gaussian.config')
