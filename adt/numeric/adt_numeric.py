from __future__ import absolute_import, unicode_literals, division, print_function
__doc__='''

This python script is developed to compute various adiabatic to diabatic transformation (ADT) quantities, like
ADT angles, ADT matrices, diabatic potential energy matrices and residue of ADT angles. It uses the python module
file, 'adt_module.so' while evaluating ADT angles and ADT matrices. User can opt any one of eight integration
paths, but the resulting diabatic potential energy matrcies are orthogonally related with each other and
therefore, dynamical observables will be independent of integration paths. For numerical calculations, user can
choose either 'HDF5' or 'text' formatted output, though the first one is fast and efficient due to low memory
requirement.

'''

__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Satyam Ravi, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''

#TODO: This code is too much scattered, reduce the complexity : C901

import os,sys
import numpy as np
from adt.numeric.adtmod import adt as fadt
from time import time
try:
    from h5py import File
    h5NotAvail  = False
except ImportError:
    h5NotAvail = True
from contextlib import contextmanager



@contextmanager
def move2dir(newdir):
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)



# This definition evaluates ADT angles, ADT matrix elements, diabatic potential energy matrix elements and residue of ADT angles

def adt_numerical(enrf, nstate, rhof, phif, path, order, outfile, logger, h5, txt, nb):
    """
        Used by the `adt` executable and shouldn't be called directly through module
    """
    start = time()
    logger.info("Reading data from files")

    # a fumction to check if the file is in numpy binary  format
    is_bin = lambda file: os.path.splitext(file)[1] in ['.npy', '.npz']
    is_h5  = lambda file: os.path.splitext(file)[1]=='.h5'

    if h5 & h5NotAvail:
        logger.error('h5py not available. Install h5py properly first.')
        sys.exit(1)
    # reading data from input files
    #NOTE: for simplicity in providing the input data using HDF5 file format,
    # input files are restricted to have just one dataset with the same name as of the file.

    if is_bin(rhof):
        rdat = np.load(rhof)
    elif is_h5(rhof):
        if h5NotAvail:
            logger.error('h5py not available. Install h5py properly first.')
            sys.exit(1)
        with File(rhof, 'r') as f:
            dset = rhof.replace('.h5', '')
            rdat = f[dset][()]
    else:
        rdat = np.loadtxt(rhof)

    if is_bin(phif):
        pdat = np.load(phif)
    elif is_h5(phif):
        with File(phif, 'r') as f:
            dset = phif.replace('.h5', '')
            pdat = f[dset][()]
    else:
        pdat = np.loadtxt(phif)
    assert rdat.shape==pdat.shape , "Mismatch in nact data"

    if enrf != None:
        if is_bin(enrf):
            enr  = np.load(enrf)
        elif is_h5(enrf):
            with File(enrf, 'r') as f:
                dset = enrf.replace('.h5', '')
                enr = f[dset][()]
        else:
            enr  = np.loadtxt(enrf)

    logger.info("Processing data")


    if nstate !=None:
        fadt.nstate = nstate
        fadt.ntau   = fadt.nstate*(fadt.nstate-1)/2
        if enrf != None: assert enr.shape[1]-2 >= fadt.nstate, \
            'Not enough data in energy file for %s states calculation.'%fadt.nstate
        assert rdat.shape[1]-2>= fadt.ntau, 'Not enough NACT data for %s states calculation.'%fadt.nstate

    else:
        fadt.ntau   = rdat.shape[1]-2
        if enrf != None:
            fadt.nstate = enr.shape[1] -2
            assert fadt.nstate*(fadt.nstate-1)/2==fadt.ntau, "Mismatch in number of states and nacts"
        else :
            i= 1
            while True:
                i+=1
                x=i*(i-1)/2
                assert x<= fadt.ntau, 'Bad number of NACTs'
                if x==fadt.ntau: break
            fadt.nstate = i
    # evaluating number of grid points, number of electronic states and number of couplings

    # set the order of the adt matrix from the provided order string if none then the default order is used
    # the order is basically a list that corresponding to the default order
    # i.e. 12,13,23 -> 1,2,3  i.e. default order
    fadt.order = getOrder(order, fadt.nstate)
    rdat   = rdat[:,:fadt.ntau+2]
    pdat   = pdat[:,:fadt.ntau+2]
    fadt.gridr  = np.unique(rdat[:,0])
    fadt.gridp  = np.unique(rdat[:,1])
    fadt.ngridr = fadt.gridr.shape[0]
    fadt.ngridp = fadt.gridp.shape[0]


    # reshaping the nonadiabatic coupling terms (NACTs) according to the number of grid points
    fadt.taur  = rdat[:,2:].reshape(fadt.ngridr, fadt.ngridp, fadt.ntau)
    fadt.taup  = pdat[:,2:].reshape(fadt.ngridr, fadt.ngridp, fadt.ntau)

    # expanding the grid points
    fadt.etaur = np.pad(fadt.taur, ((1,1),(1,1),(0,0)), "edge")
    fadt.etaup = np.pad(fadt.taup, ((1,1),(1,1),(0,0)), "edge")
    fadt.egridr= np.pad(fadt.gridr, (1,1), "reflect", reflect_type="odd")
    fadt.egridp= np.pad(fadt.gridp, (1,1), "reflect", reflect_type="odd")

    # calculation of ADT angles
    logger.info("Calculating ADT Angles on path %s"%path)
    full_angle = fadt.get_angle(fadt.ngridr, fadt.ngridp, fadt.ntau, path)
    # taking the difference between last and first value of the phi grid
    residue    = full_angle[:,-1,:] - full_angle[:,0,:] 

    full_angle = full_angle.reshape(fadt.ngridr*fadt.ngridp, fadt.ntau)
    adtAngle   = np.column_stack([rdat[:,[0,1]], full_angle])
    residue    = np.column_stack((fadt.gridr, residue))

    # calculation of ADT matrix elements
    logger.info("Calculating ADT matrix elements")
    amat =  np.apply_along_axis(fadt.amat,1,full_angle,fadt.nstate)


    if enrf != None:
        # calculation of diabatic potential energy matrix elements
        enr    = enr[:,2:fadt.nstate+2]
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
        for i in range(fadt.nstate):
            mat.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), compression="gzip")

        if enrf != None:
            logger.info("Writing Diabatic matrix elements")
            dbd = file.create_group("Diabatic Matrix elements")
            for i in range(fadt.nstate):
                dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([rdat[:,[0,1]],db[:,i,:]]), compression="gzip")

    # Writing of numerical output in '.dat' files
    if txt:
        outpath = outfile+"_%s"%path
        logger.info("Saving results in folder '%s'"%outpath)
        if not os.path.exists(outpath):os.makedirs(outpath)
        with move2dir(outpath):
            file = "Angles.dat"
            logger.info("Writing ADT Angles in '%s'"%file)
            file_write(file, adtAngle, fadt.gridr)
            # fileWrite2(file, adtAngle, fadt.ngridr, fadt.ngridp )
            logger.info("Writing ADT Angles in 'Angle_residues.dat'")
            np.savetxt('Angle_residues.dat', residue ,delimiter="\t", fmt=str("%.8f"))


            for i in range(fadt.nstate):
                file =  "Matrix_Row_%s.dat"%(i+1)
                logger.info("Writing ADT Matrix elements in '%s'"%(file))
                file_write(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), fadt.gridr )
                # fileWrite2(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]), fadt.ngridr, fadt.ngridp )

            if enrf != None:
                for i in range(fadt.nstate):
                    file =  "Diabatic_Row_%s.dat"%(i+1)
                    logger.info("Writing Diabatic Matrix elements in '%s'"%(file))
                    file_write(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]), fadt.gridr )
                    # fileWrite2(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]), fadt.ngridr, fadt.ngridp )


    # Writing of numerical output in '.npy' files
    if nb:
        outpath = outfile+"_%s"%path
        logger.info("Saving results in folder '%s'"%outpath)
        if not os.path.exists(outpath):os.makedirs(outpath)
        with move2dir(outpath):
            file = "Angles"
            logger.info("Writing ADT Angles in '%s.npy'"%file)
            np.save(file, adtAngle)
            logger.info("Writing ADT Angles in 'Angle_residues.npy'")
            np.save('Angle_residues', residue)


            for i in range(fadt.nstate):
                file =  "Matrix_Row_%s"%(i+1)
                logger.info("Writing ADT Matrix elements in '%s.npy'"%(file))
                np.save(file, np.column_stack([rdat[:,[0,1]],amat[:,i,:]]) )

            if enrf != None:
                for i in range(fadt.nstate):
                    file =  "Diabatic_Row_%s"%(i+1)
                    logger.info("Writing Diabatic Matrix elements in '%s.npy'"%(file))
                    np.save(file, np.column_stack([rdat[:,[0,1]],db[:,i,:]]) )

    logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)


# This definition is used to write the output data


def adt_numerical1d(enrf, nstate, tauf, order, outfile, logger, h5, txt, nb):
    """
        Used by the `adt` executable and shouldn't be called directly through module
    """
    start = time()
    logger.info("Reading data from files")

    # a fumction to check if the file is in numpy binary  format
    is_bin = lambda file: os.path.splitext(file)[1] in ['.npy', '.npz']
    is_h5  = lambda file: os.path.splitext(file)[1]=='.h5'

    if h5 & h5NotAvail:
        logger.error('h5py not available. Install h5py properly first.')
        sys.exit(1)
    # reading data from input files
    # NOTE: for simplicity in providing the input data using HDF5 file format,
    # input files are restricted to have just one dataset with the same name as of the file.

    if is_bin(tauf):
        taudat = np.load(tauf)
    elif is_h5(tauf):
        if h5NotAvail:
            logger.error('h5py not available. Install h5py properly first.')
            sys.exit(1)
        with File(tauf, 'r') as f:
            dset = tauf.replace('.h5', '')
            taudat = f[dset][()]
    else:
        taudat = np.loadtxt(tauf)


    if enrf != None:
        if is_bin(enrf):
            enr  = np.load(enrf)
        elif is_h5(enrf):
            with File(enrf, 'r') as f:
                dset = enrf.replace('.h5', '')
                enr = f[dset][()]
        else:
            enr  = np.loadtxt(enrf)

    logger.info("Processing data")

    # for jacobi the grid is in the form of phi, tau1 tau2 ...
    # phi is in radian
    fadt.ngrid  = taudat.shape[0]
    fadt.ntau   = taudat.shape[1]-1

    fadt.grid   = taudat[:,0]
    fadt.tau    = taudat[:,1:]

    if not np.unique(taudat[:,0]).shape == taudat[:,0].shape:
        raise Exception("Provided file doesn't have a proper 1D grid on first column")


    if nstate !=None:
        fadt.nstate = nstate
        fadt.ntau   = fadt.nstate*(fadt.nstate-1)/2
        if enrf != None: assert enr.shape[1]-1 >= fadt.nstate, \
        'Not enough data in energy file for %s states calculation.'%fadt.nstate
        assert taudat.shape[1]-1>= fadt.ntau, 'Not enough NACT data for %s states calculation.'%fadt.nstate

    else:
        fadt.ntau   = taudat.shape[1]-1
        if enrf != None:
            fadt.nstate = enr.shape[1] -1
            assert fadt.nstate*(fadt.nstate-1)/2==fadt.ntau, "Mismatch in number of states and nacts"
        else :
            i= 1
            while True:
                i+=1
                x=i*(i-1)/2
                assert x<= fadt.ntau, 'Bad number of NACTs'
                if x==fadt.ntau: break
            fadt.nstate = i
    # evaluating number of grid points, number of electronic states and number of couplings

    fadt.order = getOrder(order, fadt.nstate)


    # calculation of ADT angles
    logger.info("Calculating ADT Angles")
    full_angle = fadt.get_angle1d(fadt.ngrid,fadt.ntau)
    
    residue = full_angle[-1] - full_angle[0]
    residue = np.append(fadt.grid[-1], residue)[None]

    adtAngle   = np.column_stack([fadt.grid, full_angle])

    # calculation of ADT matrix elements
    logger.info("Calculating ADT matrix elements")
    amat = np.apply_along_axis(fadt.amat,1,full_angle,fadt.nstate)


    if enrf != None:
        # calculation of diabatic potential energy matrix elements
        enr    = enr[:,1:fadt.nstate+1]
        logger.info("Calculating Diabatic matrix elements")
        db = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat)



    # Writing of numerical output in a '.h5' file
    path = '1D'
    if h5:
        file =outfile+ "_%s.h5"%path
        logger.info("Opening HDF5 file '%s' for writing results"%file)
        file = File(file,'w')

        logger.info("Writing ADT Angles")
        ang = file.create_group("ADT Angles")
        ang.create_dataset("Angles",data=adtAngle, compression="gzip")

        ang = file.create_group("ADT Angles Residues")
        ang.create_dataset("Angles Residues",data=residue, compression="gzip")
        # logger.info('Writing ADT Angle residues')
        # ang.create_dataset("Residue", data= residue, compression="gzip")

        logger.info("Writing ADT matrix elements")
        mat = file.create_group("ADT Matrix elements")
        for i in range(fadt.nstate):
            mat.create_dataset("Row %s"%(i+1),data =np.column_stack([fadt.grid,amat[:,i,:]]), compression="gzip")

        if enrf != None:
            logger.info("Writing Diabatic matrix elements")
            dbd = file.create_group("Diabatic Matrix elements")
            for i in range(fadt.nstate):
                dbd.create_dataset("Row %s"%(i+1),data =np.column_stack([fadt.grid,db[:,i,:]]), compression="gzip")


    # Writing of numerical output in '.dat' files
    if txt:
        outpath = outfile+"_%s"%path
        logger.info("Saving results in folder '%s'"%outpath)
        if not os.path.exists(outpath):os.makedirs(outpath)
        with move2dir(outpath):
            file = "Angles.dat"
            logger.info("Writing ADT Angles in '%s'"%file)
            np.savetxt(file, adtAngle, delimiter="\t", fmt=str("%.8f"))
            logger.info("Writing ADT Angles in 'Angle_residues.dat'")
            np.savetxt('Angle_residues.dat', residue ,delimiter="\t", fmt=str("%.8f"))


            for i in range(fadt.nstate):
                file =  "Matrix_Row_%s.dat"%(i+1)
                logger.info("Writing ADT Matrix elements in '%s'"%(file))
                np.savetxt(file,  np.column_stack([fadt.grid,amat[:,i,:]]), delimiter="\t", fmt=str("%.8f"))

            if enrf != None:
                for i in range(fadt.nstate):
                    file =  "Diabatic_Row_%s.dat"%(i+1)
                    logger.info("Writing Diabatic Matrix elements in '%s'"%(file))
                    np.savetxt(file, np.column_stack([fadt.grid,db[:,i,:]]), delimiter="\t", fmt=str("%.8f"))


    # Writing of numerical output in '.npy' files
    if nb:
        outpath = outfile+"_%s"%path
        logger.info("Saving results in folder '%s'"%outpath)
        if not os.path.exists(outpath):os.makedirs(outpath)
        with move2dir(outpath):
            file = "Angles"
            logger.info("Writing ADT Angles in '%s.npy'"%file)
            np.save(file, adtAngle)
            logger.info("Writing ADT Angles in 'Angle_residues.npy'")
            np.save('Angle_residues', residue)


            for i in range(fadt.nstate):
                file =  "Matrix_Row_%s"%(i+1)
                logger.info("Writing ADT Matrix elements in '%s.npy'"%(file))
                np.save(file, np.column_stack([fadt.grid,amat[:,i,:]]) )

            if enrf != None:
                for i in range(fadt.nstate):
                    file =  "Diabatic_Row_%s"%(i+1)
                    logger.info("Writing Diabatic Matrix elements in '%s.npy'"%(file))
                    np.save(file, np.column_stack([fadt.grid,db[:,i,:]]) )

    logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)



def file_write(file, data, col):
    file = open(file, "w")
    for r in col:
        np.savetxt( file, data[data[:,0]==r] ,delimiter="\t", fmt=str("%.8f"))
        file.write("\n")
    file.flush()
    file.close()

def fileWrite2(file, data, n1, n2):
    data.shape = (n1,n2,-1)
    with open(file,'w') as f:
        for i in data:
            np.savetxt(f,i ,delimiter='\t',fmt=str('%.8f'))
            f.write('\n')


def getOrder(userOrder, state):
    if not userOrder: # no user order provided return the default
        return list(range(1,int(fadt.nstate*(fadt.nstate-1)/2+1)))

    userOrderList = userOrder.split(',')
    baseOrderList = ['{}{}'.format(i,j) for j in range(2,state+1)  for i in range(1,j)]
    check = [ordr in baseOrderList for ordr in userOrderList]
    checkR= [ordr in userOrderList for ordr in baseOrderList]
    assert all(check) and all(checkR), 'Provided order is not valid.'
    order = [baseOrderList.index(i)+1 for i in userOrderList]
    return order




# These are python modular level APIs to the ADT software package
# can be called from any python script after installation

def adt2d(grid1, grid2, nact1, nact2, energy=None, path = 1, order=None):
    '''
    Calculates ADT quantities, namely, ADT angle, residue, ADT matrix and diabatic
    matrix elements

    Parameters
    ----------
    grid1 : 1D ndarray for 1st co-ordinate grid.
    grid2 : 1D ndarray for 2nd coordinate grid.
    nact1 : 3D ndarray for component of NACT for the 1st coordinate in shape of
          `(ngrid1, ngrid2, ntau)` where ngrid1, ngrid2, ntau are the number of grid points
          for 1st coordinate , grid points for 2nd coordinate and number of NACTs respectively.
    nact2 : 3D ndarray for component of NACT for the 1st coordinate in shape same as 
          `nact1`.
    energy: 3D ndarray for energy in shape of `(ngrid1, ngrid2, nstate)` where nstate is the 
           number of energy state.


    Returns
    -------
    ADT Angle       : 3D ndarray for ADT angles in shape of `(ngrid1, ngrid2, ntau)`.
    Angle Residue   : 2D ndarray for residues of ADT angles in shape of `(ngrid1, ntau)`.
    ADT Matrix      : 4D ndarray for ADT matrix elements in shape of `(ngrid1, ngrid2, nstate, nstate)` .
    Diabatic Matrix : 4D ndarray for diabatic potential energy matrix elements 
                    in shape of `(ngrid1, ngrid2, nstate, nstate)`. Only returned when adiabatic energy is 
                    provided in the input argument.

    '''
    # quantities are given in shape of (grid1, grid2, ntau)
    fadt.gridr  = grid1
    fadt.gridp  = grid2
    fadt.ngridr = fadt.gridr.shape[0]
    fadt.ngridp = fadt.gridp.shape[0]
    fadt.ntau   = nact1.shape[2]

    # reshaping the nonadiabatic coupling terms (NACTs) according to the number of grid points
    fadt.taur  = nact1
    fadt.taup  = nact2


    if energy is None :
        i= 1
        while True:
            i+=1
            x=i*(i-1)/2
            assert x<= fadt.ntau, 'Bad number of NACTs'
            if x==fadt.ntau: break
        fadt.nstate = i
    else:
        fadt.nstate = energy.shape[2] 
        assert fadt.nstate*(fadt.nstate-1)/2==fadt.ntau, "Mismatch in number of states and nacts"

    fadt.order = getOrder(order, fadt.nstate)

    # expanding the grid points
    fadt.etaur = np.pad(fadt.taur, ((1,1),(1,1),(0,0)), "edge")
    fadt.etaup = np.pad(fadt.taup, ((1,1),(1,1),(0,0)), "edge")
    fadt.egridr= np.pad(fadt.gridr, (1,1), "reflect", reflect_type="odd")
    fadt.egridp= np.pad(fadt.gridp, (1,1), "reflect", reflect_type="odd")

    # calculation of ADT angles
    full_angle = fadt.get_angle(fadt.ngridr, fadt.ngridp, fadt.ntau, path)
    residue    = full_angle[:,-1,:] - full_angle[:,0,:]

    # calculation of ADT matrix elements

    amat =  np.apply_along_axis(fadt.amat, 1, full_angle.reshape(fadt.ngridr*fadt.ngridp, fadt.ntau) , fadt.nstate)

    if energy is None :
        amat = amat.reshape(fadt.ngridr, fadt.ngridp, fadt.nstate, fadt.nstate)
        return [full_angle, residue, amat]
    else :
        db = np.einsum("ijk,ij,ijl->ikl",amat,energy.reshape(fadt.ngridr*fadt.ngridp, fadt.nstate),amat)
        amat = amat.reshape(fadt.ngridr, fadt.ngridp, fadt.nstate, fadt.nstate)
        db = db.reshape(fadt.ngridr, fadt.ngridp, fadt.nstate, fadt.nstate)
        return [full_angle, residue, amat, db]


def adt1d(grid, taudat,energy = None, order=None,):
    '''
    Calculates ADT quantities, namely, ADT angle, residue, ADT matrix and diabatic
    matrix elements

    Parameters
    ----------
    grid   : 1D ndarray for 1st co-ordinate grid.
    nact   : 2D ndarray for NACT  in shape of
            `(ngrid, ntau)` where ngrid, ntau are the number of grid points and number of NACTs respectively.
    energy : 2D ndarray for energy in shape of `(ngrid, nstate)` where nstate is the 
            number of energy state. Optional argument.


    Returns
    -------
    ADT Angles     : 2D ndarray for ADT angles in shape of `(ngrid, ntau)`.

    ADT matrix     : 3D ndarray for ADT matrix elements in shape of `(ngrid, nstate, nstate)` .
    Angle Residues : 1D ndarray of ADT Angle residues in shape of `ntau` 
    Diabatic Matrix: 3D ndarray for diabatic potential energy matrix elements 
                    in shape of `(ngrid, nstate, nstate)`. Only returned when adiabatic energy is 
                    provided in the input argument.

    '''

    fadt.grid   = grid
    fadt.ngrid  = grid.shape[0]
    fadt.ntau   = taudat.shape[1]

    fadt.tau    = taudat

    if energy is None :
        i= 1
        while True:
            i+=1
            x=i*(i-1)/2
            assert x<= fadt.ntau, 'Bad number of NACTs'
            if x==fadt.ntau: break
        fadt.nstate = i
    else:
        fadt.nstate = energy.shape[1]
        assert fadt.nstate*(fadt.nstate-1)/2==fadt.ntau, "Mismatch in number of states and nacts"


    fadt.order = getOrder(order, fadt.nstate)

    # calculation of ADT angles
    full_angle = fadt.get_angle1d(fadt.ngrid,fadt.ntau)

    # calculation of ADT matrix elements
    residue = full_angle[-1] - full_angle[0]

    amat = np.apply_along_axis(fadt.amat,1,full_angle,fadt.nstate)


    if energy is None :
        return [full_angle, residue, amat]
    else :
        db = np.einsum("ijk,ij,ijl->ikl",amat,energy,amat)
        return [full_angle, residue, amat, db]





########################################################################################################################
