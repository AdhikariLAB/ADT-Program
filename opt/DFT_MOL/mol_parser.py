import re
import sys
import textwrap
import subprocess
import numpy as np

if sys.version_info.major>2:
    from configparser import ConfigParser as ConfigParser
else :
    from ConfigParser import SafeConfigParser as ConfigParser

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

    def runMol(self):
        try:
            runMolpro = subprocess.call(['molpro', '-n',self.optInfo['processor'], 'optg_mol.inp','--no-xml-output'])
        except:
            raise Exception('Can not run Molpro')
        if runMolpro:
            raise Exception(' Optimization failed.')

    def getResults(self):
        with open('./optg_mol.out') as f : txt = f.read()

        nAtoms =self.nAtoms
        optGeom,= re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){{6}})((?:.*\n){{{}}})'.format(nAtoms), txt)

        optGeom = [i.split()[1:] for i in filter(None, optGeom.split('\n'))]
        optGeom = np.array(optGeom, dtype=np.float)
        np.savetxt('equigeom.dat', optGeom, fmt='%.10f', delimiter='\t')

        freDat= re.findall(r' Wavenumbers \[cm-1\](.*)(?:(?:.*\n){{3}})((?:.*\n){{{}}})'.format(3*nAtoms),txt)

        freqs = np.array([j for i,_ in freDat for j in i.split() ], dtype=np.float)
        freqs = freqs[freqs!=0] # take the vibrational freqs
        np.savetxt('frequency.dat', freqs[None,:], fmt='%.10f', delimiter='\t')


        tmp= [[i.split()[1:] for i in filter(None,j.split('\n')) ] for _,j in freDat ]
        wilson = np.asarray(np.column_stack(tmp), dtype=np.float)[:,:freqs.shape[0]]
        np.savetxt('wilson.dat', wilson, fmt='%.10f', delimiter='\t')

if __name__ == "__main__":
    p = MolproOptg('mol.config')
    p.runMol()
    p.getResults()