import numpy as np
root = "/home/koushik/CODES/ADT-Program Test/TEST_RUNS/NUMERICAL_CALCULATIONS/C6H3F3+/INPUT/"
n = 6
ntau = n*(n-1)/2+1

taup = np.column_stack((np.loadtxt(root+"TAUP_1.DAT"),np.column_stack((np.loadtxt(root+"TAUP_%s.DAT"%i)[:,2:] for i in range(2,ntau)))))
taur = np.column_stack((np.loadtxt(root+"TAUR_1.DAT"), np.column_stack((np.loadtxt(root+"TAUR_%s.DAT"%i)[:,2:] for i in range(2,ntau)))))
pes =np.column_stack((np.loadtxt(root+"ADIABATICPES_1.DAT"), np.column_stack((np.loadtxt(root+"ADIABATICPES_%s.DAT"%i)[:,2:] for i in range(2,n+1)))))

g1 = np.unique(pes[:,0])

filee = open("pes1.dat","wb")
filer = open("taur1.dat","wb")
filep = open("taup1.dat","wb")
for i in g1:
    np.savetxt(filee,pes[pes[:,0]==i],delimiter="\t",fmt="%.8f")
    filee.write("\n")
    np.savetxt(filer,taur[taur[:,0]==i],delimiter="\t",fmt="%.8f")
    filer.write("\n")
    np.savetxt(filep,taup[taup[:,0]==i],delimiter="\t",fmt="%.8f")
    filep.write("\n")