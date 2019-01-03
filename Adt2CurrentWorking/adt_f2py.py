import numpy as np 
from adt_module import adt
from h5py import File
from time import time
import logging



#Create the logger
logger = logging.getLogger("ADT Numerical Program")
logger.setLevel(logging.DEBUG)
ch = logging.FileHandler("ADT.log")
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("[%(asctime)s] - %(name)s - [%(levelname)8s] - %(message)s","%Y-%m-%d %I:%M:%S %p")
ch.setFormatter(formatter)
logger.addHandler(ch)



logger.info("Starting program")
start = time()
logger.info("Reading data from files")
#%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%
rdat = np.loadtxt("INPUT/taur1.dat")
pdat = np.loadtxt("INPUT/taup1.dat")
enr  = np.loadtxt("INPUT/pes1.dat")[:,2:]
path = 6
outfile = 'ADT.h5'
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


logger.info("Processing data")
adt.gridr  = np.unique(rdat[:,0])
adt.gridp  = np.unique(rdat[:,1])
adt.ngridr = adt.gridr.shape[0]
adt.ngridp = adt.gridp.shape[0]
adt.ntau   = rdat.shape[1]-2
adt.nstate = enr.shape[1]

assert rdat.shape==pdat.shape , "Mismath in nact data"
assert adt.nstate*(adt.nstate-1)/2==adt.ntau, "Mismath in number of states and nacts"


adt.taur  = rdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)
adt.taup  = pdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)

adt.etaur = np.pad(adt.taur, ((1,1),(1,1),(0,0)), "edge")
adt.etaup = np.pad(adt.taup, ((1,1),(1,1),(0,0)), "edge")
adt.egridr= np.pad(adt.gridr, (1,1), "reflect", reflect_type="odd")
adt.egridp= np.pad(adt.gridp, (1,1), "reflect", reflect_type="odd")


#open h5 file for writing
file = File(outfile,'w')



logger.info("Calculating ADT Angles")
full_angle = adt.get_angle(adt.ngridr, adt.ngridp, adt.ntau, path).reshape(adt.ngridr*adt.ngridp, adt.ntau)
adtAngle = np.column_stack([rdat[:,[0,1]], full_angle ])


logger.info("Writing ADT Angles")
ang = file.create_group("ADT Angles")
ang.create_dataset("Angles",data=adtAngle, compression="gzip")




logger.info("Calculating ADT matrix elements")
amat =  np.apply_along_axis(adt.amat,1,full_angle,adt.nstate)

logger.info("Writing ADT matrix elements")
mat = file.create_group("ADT Matrix elements")
for i in range(adt.nstate):
    mat.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), compression="gzip")



logger.info("Calculating Diabatic matrix elements")
db = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat)

logger.info("Writing Diabatic matrix elements")
dbd = file.create_group("Diabatic Matrix elements")
for i in range(adt.nstate):
    dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],db[:,i,:]]), compression="gzip")


logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)