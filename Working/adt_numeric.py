
########################################################################################################################
#                                                                                                                      #
#    This python script is developed to compute various adiabatic to diabatic transformation (ADT) quantities, like    #
#    ADT angles, ADT matrices, diabatic potential energy matrices and residue of ADT angles. It uses the python module #
#    file, 'adt_module.so' while evaluating ADT angles and ADT matrices. User can opt any one of eight integration     #
#    paths, but the resulting diabatic potential energy matrcies are orthogonally related with each other and          #
#    therefore, dynamical observables will be independent of integration paths. For numerical calculations, user can   #
#    choose either 'HDF5' or 'text' formatted output, though the first one is fast and efficient due to low memory     #
#    requirement.                                                                                                      #
#                                                                                                                      #
#    Written by Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit     #
#    Adhikari                                                                                                          #
#                                                                                                                      #
########################################################################################################################

import os
import numpy as np 
from adt_module import adt
from h5py import File
from time import time


#  This definition evaluates ADT angles, ADT matrix elements, diabatic potential energy matrix elements and residue of ADT angles

def adt_numerical(enrf, nstate, rhof, phif, path, outfile, logger,h5, txt, nb):

    logger.info("Starting program")
    start = time()
    logger.info("Reading data from files")

    # a fumction to check if the file is in numpy binary  format
    is_bin = lambda file: os.path.splitext(file)[1] in ['.npy', '.npz']

    # reading data from input files
    if is_bin(rhof):
        rdat = np.load(rhof)
    else:
        rdat = np.loadtxt(rhof)
    if is_bin(phif):
        pdat = np.load(phif)
    else:
        pdat = np.loadtxt(phif)
    assert rdat.shape==pdat.shape , "Mismatch in nact data"

    if enrf != None:
        if is_bin(enrf):
            enr  = np.load(enrf)
        else:
            enr  = np.loadtxt(enrf)

    logger.info("Processing data")


    if nstate !=None:
        adt.nstate = nstate
        adt.ntau   = adt.nstate*(adt.nstate-1)/2
        if enrf != None: assert enr.shape[1]-2 >= adt.nstate, 'Not enough data in energy file for %s states calculation.'%adt.nstate
        assert rdat.shape[1]-2>= adt.ntau, 'Not enough NACT data for %s states calculation.'%adt.nstate

    else:
        adt.ntau   = rdat.shape[1]-2
        if enrf != None:
            adt.nstate = enr.shape[1] -2
            assert adt.nstate*(adt.nstate-1)/2==adt.ntau, "Mismatch in number of states and nacts"
        else :
        # trying automatic state calculation, not sure, checks required 
            # Using explicit solution of the quadratic equation
            # r = np.roots([1,-1,-2*adt.ntau])
            # r = r.real[abs(r.imag)<1e-6]    # 
            # adt.nstate = r[(r>1) & (r==r.astype(int))][0].astype(int)
            # Using a lazy loop to check for the number of states
            i= 1
            while True:
                i+=1
                x=i*(i-1)/2
                assert x<= adt.ntau, 'Bad number of NACTs'
                if x==adt.ntau: break
            adt.nstate = i 
            # print adt.nstate
    # evaluating number of grid points, number of electronic states and number of couplings



    rdat   = rdat[:,:adt.ntau+2]
    pdat   = pdat[:,:adt.ntau+2]
    adt.gridr  = np.unique(rdat[:,0])
    adt.gridp  = np.unique(rdat[:,1])
    adt.ngridr = adt.gridr.shape[0]
    adt.ngridp = adt.gridp.shape[0]
    

    # reshaping the nonadiabatic coupling terms (NACTs) according to the number of grid points
    adt.taur  = rdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)
    adt.taup  = pdat[:,2:].reshape(adt.ngridr, adt.ngridp, adt.ntau)

    # expanding the grid points
    adt.etaur = np.pad(adt.taur, ((1,1),(1,1),(0,0)), "edge")
    adt.etaup = np.pad(adt.taup, ((1,1),(1,1),(0,0)), "edge")
    adt.egridr= np.pad(adt.gridr, (1,1), "reflect", reflect_type="odd")
    adt.egridp= np.pad(adt.gridp, (1,1), "reflect", reflect_type="odd")

    # calculation of ADT angles
    logger.info("Calculating ADT Angles on path %s"%path)
    full_angle = adt.get_angle(adt.ngridr, adt.ngridp, adt.ntau, path)
    residue    = np.sum(full_angle, axis=1)

    full_angle = full_angle.reshape(adt.ngridr*adt.ngridp, adt.ntau)
    adtAngle   = np.column_stack([rdat[:,[0,1]], full_angle])
    residue    = np.column_stack((adt.gridr, residue))

    # calculation of ADT matrix elements
    logger.info("Calculating ADT matrix elements")
    amat =  np.apply_along_axis(adt.amat,1,full_angle,adt.nstate)


    if enrf != None:
        # calculation of diabatic potential energy matrix elements
        enr    = enr[:,2:adt.nstate+2] 
        logger.info("Calculating Diabatic matrix elements")
        db = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat)



    # Writing of numerical output in a '.h5' file
    if h5:
        file =outfile+ "_%s.h5"%path 
        logger.info("Opening HDF5 file '%s' for writing results"%file)
        file = File(file,'w')

        logger.info("Writing ADT Angles")
        ang = file.create_group("ADT Angles")
        ang.create_dataset("Angles",data=adtAngle, compression="gzip")

        logger.info('Writing ADT Angle residues')
        ang.create_dataset("Residue", data= residue, compression="gzip")

        logger.info("Writing ADT matrix elements")
        mat = file.create_group("ADT Matrix elements")
        for i in range(adt.nstate):
            mat.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), compression="gzip")

        if enrf != None:
            logger.info("Writing Diabatic matrix elements")
            dbd = file.create_group("Diabatic Matrix elements")
            for i in range(adt.nstate):
                dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],db[:,i,:]]), compression="gzip")

    # Writing of numerical output in '.dat' files
    if txt:
        outpath = outfile+"_%s"%path
        if not os.path.exists(outpath):os.makedirs(outpath)
        os.chdir(outpath)
        file = "Angles.dat"
        logger.info("Writing ADT Angles in '%s'"%file)
        file_write(file, adtAngle, adt.gridr)
        logger.info("Writing ADT Angles in 'Angle_residues.dat'")
        np.savetxt('Angle_residues.dat', residue ,delimiter="\t", fmt="%.8f")

        
        for i in range(adt.nstate):
            file =  "Matrix_Row_%s.dat"%(i+1)
            logger.info("Writing ADT Matrix elements in '%s'"%(file))
            file_write(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), adt.gridr )

        if enrf != None:
            for i in range(adt.nstate):
                file =  "Diabatic_Row_%s.dat"%(i+1)
                logger.info("Writing Diabatic Matrix elements in '%s'"%(file))
                file_write(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]), adt.gridr )


    # Writing of numerical output in '.npy' files
    if nb:
        outpath = outfile+"_%s"%path
        if not os.path.exists(outpath):os.makedirs(outpath)
        os.chdir(outpath)
        file = "Angles"
        logger.info("Writing ADT Angles in '%s.npy'"%file)
        np.save(file, adtAngle)
        logger.info("Writing ADT Angles in 'Angle_residues.npy'")
        np.save('Angle_residues', residue)

        
        for i in range(adt.nstate):
            file =  "Matrix_Row_%s"%(i+1)
            logger.info("Writing ADT Matrix elements in '%s.npy'"%(file))
            np.save(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]) )

        if enrf != None:
            for i in range(adt.nstate):
                file =  "Diabatic_Row_%s"%(i+1)
                logger.info("Writing Diabatic Matrix elements in '%s.npy'"%(file))
                np.save(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]) )



    logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)


# This definition is used to write the output data
 
def file_write(file, data, col):
    file = open(file, "w")
    for r in col:
        np.savetxt( file, data[data[:,0]==r] ,delimiter="\t", fmt="%.8f")
        file.write("\n")

########################################################################################################################
