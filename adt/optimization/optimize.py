import re
import sys
import textwrap
import subprocess
import numpy as np

if sys.version_info.major>2:
    from configparser import ConfigParser as ConfigParser
else :
    from ConfigParser import SafeConfigParser as ConfigParser

class Gaussianoptg():

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
    
    def runGauss(self):
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
        np.savetxt('equigeom.dat', optGeom, fmt='%.10f', delimiter='\t')
        
        freDat= re.findall(r' Frequencies --(.*)(?:(?:.*\n){{5}})((?:.*\n){{{}}})'.format(atomNo),txt)
        freq = np.array( [i.split() for i,_ in freDat], dtype=np.float64).reshape(-1)
        np.savetxt('frequency.dat', freq[None,:], fmt='%.10f', delimiter='\t')
        
        wilson = np.array([[i.split()  for i in filter(None,j.split('\n')) ] for _,j in freDat ], dtype=np.float64)
        # remove columns regarding the atom number and atom type (one row for one atom)
        wilson = np.column_stack(wilson[...,2:])
        # structure the array, as required by the ADT code
        res = np.vstack(np.transpose(wilson.reshape(atomNo,-1,3),(0,2,1)))
        np.savetxt('wilson.dat', res, fmt='%.10f', delimiter='\t')

if __name__ == "__main__":
    p = Gaussianoptg('gaussian.config')
