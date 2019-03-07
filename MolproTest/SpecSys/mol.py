import subprocess,os,shutil
import numpy as np

# build a template

# equilibrium step

scrdir = '/tmp/koushik'

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
tmp= " {}\n".format(natoms)
tmp+= "This geometry is for Equilibirium geometry.\n"
for i in range(natoms):
    tmp += "{},{},{},{}\n".format(natomsN[i], *tmpGeom[i])


with open("geom.xyz","w") as f:
    f.write(tmp)


exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','equi.com'])
if exitcode==0:
    print "success"

with open('molequienr.res',"r") as f:
    dat = f.read().strip().split("\n")[1:]
    print "{}\t{}".format(*dat)



# create a backup of the current wavefunction


# os.remove('equi.out')
# os.remove('equi.xml')


rho_grid = np.arange(0.1,5.1, 0.1)
phi_grid = np.arange(1,181, 2)


txte= ''
txtn =''




for phi in phi_grid[:1]:
    for rho in rho_grid[:11]:

        shutil.copy('molpro.wfu',scrdir)

        Q      = np.array([ rho*np.cos(np.deg2rad(phi)), rho*np.sin(np.deg2rad(phi)), 0])

        dsGeom = np.dot(wilson, np.sqrt(hbar/freq)*Q)/np.sqrt(mass)

        curGeom = equiGeom+dsGeom

        curGeom = curGeom.reshape(natoms, 3)
        tmp= " {}\n".format(natoms)
        tmp+= "This geometry is for rho = {}, phi={}\n".format(rho,phi)
        for i in range(natoms):
            tmp += "{},{},{},{}\n".format(natomsN[i], *curGeom[i])


        with open("geom.xyz","w") as f:
            f.write(tmp)


        exitcode = subprocess.call(['molpro',"-d", scrdir,'-W .','grid.com'])
        if exitcode==1:
            continue

        # if job is successful, then update the bcakup
        # shutil.copy(scrdir+'/molpro.wfu', scrdir+'/backup.wfu')

        with open('molenr.res',"r") as f:
            dat = f.read().strip().split("\n")[1:]
            print 
            txte += "{:.2f}\t{}{}{}\n".format(rho, phi, *dat)
        with open('nac1.res',"r") as f:
            dat = f.read().replace("D","E").strip().split("\n")[1:]

            dat = np.array([float(j) for i in dat for j in i.strip().split() ])



            mtau  = dat*angtobohr*np.sqrt(hbar/mass)
            tauq1 = np.sum(mtau*wilson[:,0]/np.sqrt(freq[0]))
            tauq2 = np.sum(mtau*wilson[:,1]/np.sqrt(freq[1]))

            taurho = tauq1*np.cos(np.deg2rad(phi)) + tauq2*np.sin(np.deg2rad(phi))
            tauphi = -rho*tauq1*np.sin(np.deg2rad(phi)) + rho*np.cos(np.deg2rad(phi))*tauq2 
            print rho , phi, taurho, tauphi 


