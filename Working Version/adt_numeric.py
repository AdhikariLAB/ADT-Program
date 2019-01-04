import numpy as np 
from adt_module import adt
from h5py import File
from time import time


def adt_numerical(enrf, rhof, phif, path, outfile, logger,h5, txt):

    logger.info("Starting program")
    start = time()
    logger.info("Reading data from files")

    rdat = np.loadtxt(rhof)
    pdat = np.loadtxt(phif)
    enr  = np.loadtxt(enrf)[:,2:]
    path = path

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






    logger.info("Calculating ADT Angles on path %s"%path)
    full_angle = adt.get_angle(adt.ngridr, adt.ngridp, adt.ntau, path).reshape(adt.ngridr*adt.ngridp, adt.ntau)
    adtAngle = np.column_stack([rdat[:,[0,1]], full_angle ])


    logger.info("Calculating ADT matrix elements")
    amat =  np.apply_along_axis(adt.amat,1,full_angle,adt.nstate)



    logger.info("Calculating Diabatic matrix elements")
    db = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat)




    if h5:
        file =outfile+ "_%s.h5"%path 
        logger.info("Opening HDF5 file '%s' for writing results"%file)
        file = File(file,'w')

        logger.info("Writing ADT Angles")
        ang = file.create_group("ADT Angles")
        ang.create_dataset("Angles",data=adtAngle, compression="gzip")

        logger.info("Writing ADT matrix elements")
        mat = file.create_group("ADT Matrix elements")
        for i in range(adt.nstate):
            mat.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), compression="gzip")

        logger.info("Writing Diabatic matrix elements")
        dbd = file.create_group("Diabatic Matrix elements")
        for i in range(adt.nstate):
            dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],db[:,i,:]]), compression="gzip")

    if txt:
        file = outfile+ "_Angles.dat"
        logger.info("Writing ADT Angles in '%s'"%file)
        file_write(file, adtAngle, adt.gridr)

        
        for i in range(adt.nstate):
            file = outfile + "_Matrix_Row_%s.dat"%(i+1)
            logger.info("Writing ADT Matrix elements in '%s'"%(file))
            file_write(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), adt.gridr )


        for i in range(adt.nstate):
            file = outfile + "_Diabatic_Row_%s.dat"%(i+1)
            logger.info("Writing Diabatic Matrix elements in '%s'"%(file))
            file_write(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]), adt.gridr )


    logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)


def file_write(file, data, col):
    file = open(file, "w")
    for r in col:
        np.savetxt( file, data[data[:,0]==r] ,delimiter="\t", fmt="%.8f")
        file.write("\n")