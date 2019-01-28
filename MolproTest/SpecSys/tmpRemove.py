
import numpy as np



natoms = 3
natomsN = ["N","O","O"]

angtobohr = 1.88973
hbar = 0.063508
wilson = np.loadtxt("/home/koushik/CODES/ADT-Program_Test/MolproTest/OPT/wilson.dat")
freq = np.loadtxt("/home/koushik/CODES/ADT-Program_Test/MolproTest/OPT/freq.dat")
freq = freq*0.001883651


equiGeom = np.loadtxt("/home/koushik/CODES/ADT-Program_Test/MolproTest/OPT/equiGeom.dat")
equiGeom = equiGeom.reshape(9,1)[:,0]


mass =  np.array([ 14.006700,14.006700,14.006700, 15.999400 ,15.999400 ,15.999400 ,15.999400 ,15.999400 ,15.999400 ])




########## Equilibrium geometry step



tmpGeom = equiGeom.reshape(natoms,3)

# create a backup of the current wavefunction



rho_grid = np.arange(0.1,5.1, 0.1)
phi_grid = np.arange(1,181, 2)



for phi in phi_grid[:1]:
    for rho in rho_grid[:1]:


        Q      = np.array([ rho*np.cos(np.deg2rad(phi)), rho*np.sin(np.deg2rad(phi)), 0])

        dsGeom = np.dot(wilson, np.sqrt(hbar/freq)*Q)/np.sqrt(mass)
        print dsGeom.reshape(natoms,3)


nModes = 3
wilson = wilson.reshape(natoms, 3,  nModes)
mass =  np.array([ 14.006700, 15.999400 ,15.999400 ])
for phi in phi_grid[:1]:
    for rho in rho_grid[:1]:


        Q      = np.array([ rho*np.cos(np.deg2rad(phi)), rho*np.sin(np.deg2rad(phi)), 0])

        dsGeom = np.einsum('ijk,k,k,i->ij',wilson,np.sqrt(hbar/freq),Q,np.sqrt(1/mass))


        print dsGeom
