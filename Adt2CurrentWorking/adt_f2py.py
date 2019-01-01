import numpy as np 
from adt_module import adtfuncs as adt
import h5py
from time import time
start = time()
rdat = np.loadtxt("../INPUT/taur1.dat")
pdat = np.loadtxt("../INPUT/taup1.dat")
enr  = np.loadtxt("../INPUT/pes1.dat")[:,2:]



adt.gridr  = np.unique(rdat[:,0])
adt.gridp  = np.unique(rdat[:,1])
adt.ngridr = adt.gridr.shape[0]
adt.ngridp = adt.gridp.shape[0]
adt.ntau   = rdat.shape[1]-2
adt.nstate = enr.shape[1]

adt.taur  = rdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)
adt.taup  = pdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)

adt.etaur = np.pad(adt.taur, ((1,1),(1,1),(0,0)), "edge")
adt.etaup = np.pad(adt.taup, ((1,1),(1,1),(0,0)), "edge")
adt.egridr= np.pad(adt.gridr, (1,1), "reflect", reflect_type="odd")
adt.egridp= np.pad(adt.gridp, (1,1), "reflect", reflect_type="odd")
path = 6

full_angle = adt.get_angle(adt.ngridr, adt.ngridp, adt.ntau, path).reshape(adt.ngridr*adt.ngridp, adt.ntau)
adtAngle = np.column_stack([rdat[:,[0,1]], full_angle ])



file = h5py.File('ADT.h5','w')
ang = file.create_group("ADT Angles")
ang.create_dataset("Angles",data=adtAngle, compression="gzip")


amat =  np.apply_along_axis(adt.amat,1,full_angle,adt.nstate)

mat = file.create_group("ADT Matrix elements")
for i in range(adt.nstate):
    mat.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), compression="gzip")


#calculate diabatic matrices
db = np.stack(np.einsum("ijk,ij,ijl->ikl",amat,enr,amat))
dbd = file.create_group("Diabatic Matrix elements")
for i in range(adt.nstate):
    dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],db[:,i,:]]), compression="gzip")


print "The program took %s seconds"time()-start