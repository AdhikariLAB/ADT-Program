import re 
import numpy as np 

with open('./optg_mol.out') as f : txt = f.read()

nAtoms =self.nAtoms# int(re.findall(r' END OF GEOMETRY OPTIMIZATION.(?:(?:.*\n){4})\s+(\d+)', txt)[0])
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
