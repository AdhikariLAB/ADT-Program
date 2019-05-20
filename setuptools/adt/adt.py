
__doc__='''

This python file parses the command line arguments associated with 'adt' command. It uses subparser either to
devise analytic functional forms of several adiabatic to diabatic transformation (ADT) quantities or to solve the
stiff ADT equations for any `N' coupled electronic states. While carrying out symbolic manipulation, it employes
the definitions of adt_analytic.py. On the other hand, adt_numeric.py is involved for numerical calculation of
ADT angles, ADT matrcies, diabatic potential energy matrices and residue of ADT angles. In order to monitor the
progress of a job, a auto-generated log file, 'ADT.log' is created during the execution.

Any user can easily get the help message by typing 'adt -h' for the overall outline of this program. On the other
hand, 'adt ana -h' or 'adt num -h' can be executed for more specific informations about analytical or numerical
jobs.

'''
__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''


import os
import sys
import logging
import textwrap
import argparse
from analytic.adt_analytic import adt_analytical
from molpro.adt_molpro import Scattering, Spectroscopic, Jacobi
import ConfigParser

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('\033[91mError: %s\n\033[0m' % message)
        self.print_help()
        sys.exit(2)



def make_logger(log_name):
    #Create the logger
    logger = logging.getLogger(log_name)
    logger.setLevel(logging.DEBUG)
    # file = os.path.join(os.path.expanduser('~'), 'ADT.log' )
    fh = logging.FileHandler("ADT.log")
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] - %(name)22s - [%(levelname)6s] - %(message)s","%d-%m-%Y %I:%M:%S %p")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def main():
    #main parser
    parser = CustomParser(
        prog="adt",
        formatter_class=argparse.RawTextHelpFormatter,
        description = textwrap.dedent('''
    A generalised ADT program for analytical and numerical calculation. This is applicable for any
    'N' electronic state sub-Hilbert space. The analytical segment can be used to generate symbolic
    expressions of eigth adiabatic to diabatic transformation (ADT) quantities (elements of adiabatic
    potential energy matrix, elements of nonadiabatic coupling matrix, ADT matrix elements, partially
    substituted ADT equations, completely substituted ADT equations, elements of coefficient matrix of
    gradient of ADT angles, elements  of coefficient matrix of nonadiabatic coupling terms (NACTs) and
    the diabatic potential energy matrix elements) for any arbitrary number of coupled electronic
    states. On the other hand, the numerical portion computes ADT angles, ADT matrix elements, diabatic
    potential energy matrix elements and residue of ADT angles for any 'N' coupled electronic states
    with multiple degrees of freedom. Any user can solve the differential equations along eight
    different paths over the nuclear configuration space (CS).''')
        )

    #adding subparsers
    subparsers = parser.add_subparsers(title='Available sub-commands', dest="choice", help="choose from one of these")


    #subparser for analytical jobs
    analytical = subparsers.add_parser("ana",
        formatter_class=argparse.RawTextHelpFormatter,
        description = textwrap.dedent('''
    Devise analytical expressions of any one of the ADT quantities for a given number of states'''),
        help ="Formulate analytical expressions")
    analytical_required = analytical.add_argument_group("Required arguments")



    #subparser for numerical jobs
    numeric = subparsers.add_parser("num",
        formatter_class=argparse.RawTextHelpFormatter,
        description = textwrap.dedent('''
    Calculate ADT angle and diabatic potential energy matrix for a given number of electronic states along a specific path.
    '''),
        help= "Calculate ADT angle and diabatic potential energy matrix")
    numeric_required = numeric.add_argument_group("Required arguments")


    molpro = subparsers.add_parser('mol',
        formatter_class=argparse.RawTextHelpFormatter,
        description = textwrap.dedent('''
    Calculate the Adiabatic potential energy surfaces(PESs) and noadiabatic coupling matrix(NACM)
    and subsequently calculate the Numerical quantities using the numerical section.'''),
        help= "Run MOLPRO and calculate ADT angles and diabatic surfaces and couplings")
    molpro_required = molpro.add_argument_group("Required arguments")


    #adding options for analytical jobs
    analytical_required.add_argument("-nstate",
                                    type=int,
                                    help="Number of states", required=True)
    analytical         .add_argument("-anajob",
                                    type=int,
                                    help="Specify the type of expression (default: %(default)s - completely substituted ADT equation) ",
                                    choices=range(1,9),
                                    metavar="[1-8]",
                                    default=5)



    #adding options for numerical jobs
    numeric_required.add_argument("-nfile1",
                        type     = str,
                        help     = "Specify the NACT file along first coordinate. \nThis one NACT file is sufficient for evaluating ADT for a 1D grid of geometries. \nBut for calulating a 2D ADT user must also provide another NACT file \n(usinf '-nfile2') along the other coordinate.",
                        metavar  = "FILE",
                        required = True)
    numeric.add_argument("-nfile2",
                        type     = str,
                        help     = "Specify the NACT file along second coordinate. \nRequired for calculating ADT over a 2D grid of geometries.",
                        metavar  = "FILE")
    numeric.add_argument("-intpath",
                        type    = int,
                        help    = "Specify the path for calculation (default: %(default)s).\n ",
                        choices = range(1,9),
                        metavar = "[1-8]",
                        default = 1)
    numeric.add_argument("-efile",
                        type    = str,
                        help    = "Specify the Energy file for calculating the diabatic potential energy matrix elements.\n ",
                        metavar = "FILE")
    numeric.add_argument('-nstate',
                        type = int,
                        help = "Specify the number of states to do the calculation.\nBy default it includes all the data for calculation.\n  ")
    numeric.add_argument("-ofile",
                        type    = str,
                        help    = "Specify the output folder/file name (w/o extension) (default: %(default)s).\n ",
                        metavar = "FILE",
                        default="'ADT_numeric'")
    numeric.add_argument("-n",
                        type    =str,
                        help="Specify number of OpenMP threads to use for parallel calculation. \nApplicable only when installed using OpenMP support.\n(default: Maximum avilable threads)\n ",
                        default= False)
    numeric.add_argument("-h5",
                        action = 'store_true',
                        help   = "Write results in a HDF5 file (.h5). \nFast IO, smaller file size and hierarchical filesystem-like data format,\npreferable for saving and sharing large datasets in an organised way.\n " )
    numeric.add_argument("-nb",
                        action = 'store_true',
                        help   = 'Write results in Numpy binary file (.npy). \nPreferable when working with numpy for its much faster IO and easy portability.\n ')
    numeric.add_argument("-txt" ,
                        action = "store_true",
                         help="Write results in a text file.(default behaviour).")
    numeric.set_defaults(h5=False,txt=False, nb = False)



    molpro_required.add_argument('-config',
                        type    = str,
                        metavar = "FILE",
                        required = True,
                        help='Specify the molpro configuration file containing\nthe necessary keywords of MOLPRO software.\n  ')
    molpro_required.add_argument('-atomfile',
                        type    = str,
                        metavar = "FILE",
                        required = True,
                        help='Specify the information file constituting atomic\nsymbols and atomic masses. \n ')
    molpro.add_argument('-geomfile',
                        type    = str,
                        metavar = "FILE",
                        default='geomfile.dat',
                        help='Specify the geometry file containing the initial\ngrid point in "xyz" format (in Angstrom)\n(default: %(default)s). \n(Required for spectroscopic system). \n ')
    molpro.add_argument('-freqfile',
                        type    = str,
                        metavar = "FILE",
                        default = 'frequency.dat',
                        help='Specify the frequency information file, where\nfrequencies of normal modes are written in cm-1\n(default: %(default)s). \n(Required for spectroscopic system). \n ')
    molpro.add_argument('-wilsonfile',
                        type    = str,
                        metavar = "FILE",
                        default = 'wilson.dat',
                        help='Specify the filename containing the Wilson matrix\nof a molecular species (default: wilson.dat).\n(default: %(default)s).\n(Required for spectroscopic system). \n ')
    molpro.add_argument("-intpath",
                        type    = int,
                        help    = "Specify the path for calculation (default: %(default)s).\n ",
                        choices = range(1,9),
                        metavar = "[1-8]",
                        default = 1)
    molpro.add_argument("-ofile",
                        type    = str,
                        help    = "Specify the output file name (w/o extension) (default: %(default)s).\n ",
                        metavar = "FILE",
                        default="'ADT_numeric'")
    molpro.add_argument("-n",
                         type=str,
                         help="Specify number of OpenMP threads to use for parallel calculation. \nApplicable only when installed using OpenMP support.\n(default: Maximum avilable threads)\n ",
                         default=False)
    molpro.add_argument("-h5",
                        action = 'store_true',
                        help   = "Write results in a HDF5 file (.h5).\nFast IO, smaller file size and hierarchical filesystem-like data format,\npreferable for saving and sharing large datasets in an organised way.\n " )
    molpro.add_argument("-nb",
                        action = 'store_true',
                        help   = 'Write results in Numpy binary file (.npy). \nPreferable when working with numpy for its much faster IO and easy portability.\n ')
    molpro.add_argument("-txt" ,
                        action = "store_true",
                        help="Write results in a text file. (default behaviour).")
    molpro.set_defaults(h5=False,txt=False, nb = False)


    args = parser.parse_args()


    #generation of analytic expressions according to command line arguments
    if args.choice =="ana":
        state  = args.nstate
        path   = args.anajob
        logger = make_logger("ADT Analytical program")

        logger.info('''Starting Analytical program.

        Number of states   : {}
        Integration path   : {}
        '''.format(state, path))
        try:
            adt_analytical(state, path, logger)
            print("Log saved in 'ADT.log'.")
        except Exception as e:
            logger.error("Program failed. %s\n"%e+"-"*121)
            print("Program failed. %s"%e)



    #calculation of ADT quantities according to command line arguments
    if args.choice == "num":
        path    = args.intpath
        enrf    = args.efile
        nstate  = args.nstate
        rhof    = args.nfile1
        phif    = args.nfile2
        outfile = args.ofile.strip("'")
        h5      = args.h5
        txt     = args.txt
        nb = args.nb
        threads = args.n

        if (h5==False and txt== False and nb==False ) : txt=True

        ffrmt = []
        if h5: ffrmt.append('HDF5')
        if txt: ffrmt.append('Plain txt')
        if nb: ffrmt.append('Numpy Binary format')
        ffrmt = ', '.join(ffrmt)

        if not enrf : nstatet = None
        if enrf and not nstate : nstatet = 'All'
        else :nstatet = nstate
        
        # adt1D true means an 1D adt will be done
        adt2D = True if phif else False



        if adt2D :
            logger = make_logger("ADT Numerical Program")
            tmpLog = '''Starting Numerical program for 2D ADT. 

            Energy File          : {}
            NACT 1 File          : {}
            NACT 2 File          : {}
            Number of states     : {}
            Integration path     : {}
            Output file/folder   : {}
            Output file format   : {}
            '''.format(enrf, rhof, phif, nstatet, path, outfile, ffrmt)
            if threads:
                tmpLog += 'OpenMP threads       : {}'.format(threads)
                # set opnemp environment variable to spwan threads
                os.environ['OMP_NUM_THREADS'] = threads
            logger.info(tmpLog)
            
            # importing it here, just so that `OMP_NUM_THREADS` can take effect
            from numeric.adt_numeric import adt_numerical

            try:
                adt_numerical(enrf, nstate, rhof, phif, path, outfile, logger, h5, txt, nb)
                print("Log saved in 'ADT.log'.")
            except Exception as e:
                logger.error("Program failed. %s\n"%e+"-"*121)
                print("Program failed. %s"%e)
        else :
            logger = make_logger("ADT Numerical Program")
            tmpLog = '''Starting Numerical program for 1D ADT. 

            Energy File          : {}
            NACT   File          : {}
            Number of states     : {}
            Output file/folder   : {}
            Output file format   : {}
            '''.format(enrf, rhof, nstatet, outfile, ffrmt)

            #!!! no threaded part for 1D adt
            # if threads:
            #     tmpLog += 'OpenMP threads       : {}'.format(threads)
            #     # set opnemp environment variable to spwan threads
            #     os.environ['OMP_NUM_THREADS'] = threads
            logger.info(tmpLog)
            
            # importing it here, just so that `OMP_NUM_THREADS` can take effect
            from numeric.adt_numeric import adt_numerical1d

            try:
                adt_numerical1d(enrf, nstate, rhof, outfile, logger, h5, txt, nb)
                print("Log saved in 'ADT.log'.")
            except Exception as e:
                logger.error("Program failed. %s\n"%e+"-"*121)
                print("Program failed. %s"%e)




    if args.choice == 'mol':
        # arguments related to running molpro 
        configfile = args.config
        atomfile   = args.atomfile
        geomfile   = args.geomfile
        freqfile   = args.freqfile
        wilsonfile = args.wilsonfile


        # arguments related to numerical calculation
        path       = args.intpath
        outfile = args.ofile.strip("'")
        h5      = args.h5
        txt     = args.txt
        nb = args.nb
        
        if (h5 == False and txt == False and nb == False): txt = True

        logger = make_logger("ADT Molpro Program")

        ffrmt = []
        if h5: ffrmt.append('HDF5')
        if txt: ffrmt.append('Plain txt')
        if nb: ffrmt.append('Numpy Binary format')
        ffrmt = ', '.join(ffrmt)



        scf = ConfigParser.SafeConfigParser()
        scf.read(configfile)
        sysType = dict(scf.items('sysInfo'))['type']


        infoText = "Starting molpro jobs. "
        outfoText='''Molpro jobs completed. Data saved in following files:
                Energy File : energy_mod.dat'''

        try:
            if sysType == 'spectroscopic':
                infoText +='''
                System type               : Spectroscopic
                Co-ordinate type          : Normal Modes
                Molpro Config file        : {}
                Atom Info file            : {}
                Equilibrium Geometry file : {}
                Frequency Info file       : {}
                Wilson Matrix file        : {}

                '''.format(configfile, atomfile, geomfile, freqfile, wilsonfile)
                logger.info(infoText + "\nCheck 'adt_molpro.log' for progress...")
                jobRunner = Spectroscopic(scf, atomfile, geomFile, freqfile, wilsonFile)
                trFile = 'tau_rho_mod.dat'
                outfoText += '''
                NACT 1 File : tau_rho_mod.dat
                NACT 2 File : '''



            elif sysType == 'scattering_hyper':
                infoText +='''
                System type               : Scattering
                Co-ordinate type          : Hyperspherical
                Molpro Config file        : {}
                Atom Info file            : {}

                '''.format(configfile, atomfile)
                logger.info(infoText + "\nCheck 'adt_molpro.log' for progress...")
                jobRunner = Scattering(scf, atomfile)
                trFile = 'tau_rho_mod.dat'
                outfoText += '''
                NACT 1 File : tau_rho_mod.dat
                NACT 2 File : '''

            elif sysType == 'scattering_jacobi':
                infoText +='''
                System type               : Scattering
                Co-ordinate type          : Jacobi
                Molpro Config file        : {}
                Atom Info file            : {}

                '''.format(configfile, atomfile)

                logger.info(infoText + "\nCheck 'adt_molpro.log' for progress...")
                jobRunner = Jacobi(scf, atomfile)
                outfoText += '''
                NACT File : '''

            else :
                raise Exception('Not a proper system type')
            jobRunner.runMolpro()
            logger.info(outfoText+'tau_phi_mod.dat')

        except Exception as e:
            logger.error("Program failed in molpro job. %s\n"%e+"-"*121)
            sys.exit(1)

        try:
            from numeric.adt_numeric import adt_numerical, adt_numerical1d
            if sysType in ['spectroscopic' , 'scattering_hyper']:
                logger.info('''Starting Numerical calculation
                Integration Path   : {}
                Output file/folder : {}
                Output file format : {}
                '''.format(path, outfile, ffrmt))
                adt_numerical('energy_mod.dat', None, trFile, 'tau_phi_mod.dat', path, outfile, logger, h5, txt, nb)
            else :
                logger.info('''Starting Numerical calculation
                Output file/folder : {}
                Output file format : {}
                '''.format(path, outfile, ffrmt))
                adt_numerical1d('energy_mod.dat', None, 'tau_phi_mod.dat', outfile, logger, h5, txt, nb)
        except Exception as e:
            logger.error("Program failed in numerical calculation. %s\n"%e+"-"*121)

    #######################################################################################################################

if __name__ == "__main__":
    main()
