
########################################################################################################################
#                                                                                                                      #
#    This python file parses the command line arguments associated with 'adt' command. It uses subparser either to     #
#    devise analytic functional forms of several adiabatic to diabatic transformation (ADT) quantities or to solve the #
#    stiff ADT equations for any `N' coupled electronic states. While carrying out symbolic manipulation, it employes  #
#    the definitions of adt_analytic.py. On the other hand, adt_numeric.py is involved for numerical calculation of    #
#    ADT angles, ADT matrcies, diabatic potential energy matrices and residue of ADT angles. In order to monitor the   #
#    progress of a job, a auto-generated log file, 'ADT.log' is created during the execution.                          #
#                                                                                                                      #  
#    Any user can easily get the help message by typing 'adt -h' for the overall outline of this program. On the other #
#    hand, 'adt ana -h' or 'adt num -h' can be executed for more specific informations about analytical or numerical   #
#    jobs.                                                                                                             #
#                                                                                                                      #
#    Written by Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and              #
#    Satrajit Adhikari                                                                                                 #
#                                                                                                                      #
########################################################################################################################

import argparse,sys,logging
from adt_numeric import *
from adt_analytic import *
import textwrap



#main parser
parser = argparse.ArgumentParser(
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
Calculate numerical results for a given number of electronic states along a specific path'''),
    help= "Calculate numerical results")
numeric_required = numeric.add_argument_group("Required arguments")



#adding options for analytical jobs
analytical_required.add_argument("-nstate", type=int, help="Number of states", required=True)
analytical         .add_argument("-anajob", type=int, help="Specify the type of expression (default: %(default)s - completely substituted ADT equation) ", choices=range(1,9), metavar="[1-8]", default=5)



#adding options for numerical jobs
numeric         .add_argument("-intpath", type=int, help="Specify the path for calculation (default: %(default)s).\n ", choices=range(1,9),metavar="[1-8]", default=1)
numeric         .add_argument("-efile",   type=str, help="Specify the Energy file for calculating the diabatic potential energy matrix elements.\n ", metavar="FILE")
numeric         .add_argument('-nstate',  type=int, help="Specify the number of states to do the calculation.\nBy default it includes all the data for calculation.\n  ")
numeric         .add_argument("-ofile",   type=str, help="Specify the output file name (w/o extension) (default: %(default)s).\n ", metavar="FILE", default="'ADT_numeric'")
numeric_required.add_argument("-nfile1",  type=str, help="Specify the NACT file along first coordinate.\n ", metavar="FILE", required=True)
numeric_required.add_argument("-nfile2",  type=str, help="Specify the NACT file along second coordinate.", metavar="FILE", required=True)
numeric.add_argument("-h5", action='store_true',    help="Write results in a HDF5 file (.h5). (default behaviour).\nFast IO, smaller file size and hierarchical filesystem-like data format, \npreferable for saving and sharing large datasets in an organised way.\n " )
numeric.add_argument("-nb", action='store_true',    help='Write results in Numpy binary file (.npy). \nPreferable when working with numpy for its much faster IO and easy portability.\n ')
numeric.add_argument("-txt" ,action="store_true",   help="Write results in a text file.")
numeric.set_defaults(h5=False,txt=False, nb = False)




def make_logger(log_name):
    #Create the logger
    logger = logging.getLogger(log_name)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler("ADT.log")
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(asctime)s] - %(name)22s - [%(levelname)6s] - %(message)s","%Y-%m-%d %I:%M:%S %p")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


#collecting arguments
args = parser.parse_args()


#generation of analytic expressions according to command line arguments
if args.choice =="ana":
    state = args.nstate 
    path  = args.anajob
    logger = make_logger("ADT Analytical program")
    try:
        adt_analytical(state, path, logger)
        print("Log saved in 'ADT.log'.")
    except Exception as e:
        logger.error("Program failed\n"+"-"*121)
        print("Program failed. %s"%e)



#calculation of ADT quantities according to command line arguments
if args.choice == "num":
    path = args.intpath 
    enrf = args.efile 
    nstate = args.nstate
    rhof = args.nfile1 
    phif = args.nfile2
    outfile = args.ofile.strip("'")
    h5 = args.h5
    txt = args.txt
    nb = args.nb
    if (h5==False and txt== False and nb==False ) : h5=True 

    logger = make_logger("ADT Numerical Program")
    try:
        adt_numerical(enrf,nstate, rhof, phif, path, outfile, logger, h5, txt, nb)
        print("Log saved in 'ADT.log'.")
    except Exception as e:
        logger.error("Program failed\n"+"-"*121)
        print("Program failed. %s"%e)

#######################################################################################################################
