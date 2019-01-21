
############################################################################################################################
# Authors are Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
############################################################################################################################

import argparse,sys,logging
from adt_numeric import *
from adt_analytic import *



#main parser
parser = argparse.ArgumentParser(
    prog="ADT",
    formatter_class=argparse.RawTextHelpFormatter,
    description = "A generalised ADT program for analytical and numerical calculation."
    )

#adding subparsers
subparsers = parser.add_subparsers(title='Available sub-commands', dest="choice", help="choose from one of these")


#subparser for analytical jobs
analytical = subparsers.add_parser("ana",
    formatter_class=argparse.RawTextHelpFormatter,
    description = "Calculate different analytical expression for a given number of states",
    help ="Calculate different analytical expression")
analytical_required = analytical.add_argument_group("Required arguments")



#subparser for numerical jobs
numeric = subparsers.add_parser("num", 
    formatter_class=argparse.RawTextHelpFormatter,
    description = "Calculate numerical results for a given number of states",
    help= "Calculate numerical results")
numeric_required = numeric.add_argument_group("Required arguments")



#adding options for analytical jobs
analytical_required.add_argument("-n", type=int, help="Number of states", required=True)
analytical         .add_argument("-p", type=int, help="Specify the type of expression (default: %(default)s - Full ADT equation) ", choices=range(1,9), metavar="[1-8]", default=5)



#adding options for numerical jobs
numeric         .add_argument("-p", type=int, help="Specify the path for calculation (default: %(default)s)", choices=range(1,9),metavar="[1-8]", default=1)
numeric         .add_argument("-o", type=str, help="Specify the output file name (w/o extension) (default: %(default)s)", metavar="FILE", default="'ADT_numeric'")
# numeric         .add_argument("-h5", default=True, metavar="True/False", help="Write results in a HDF5 file (fast and efficient)  (default: %(default)s)",action='store_false')
numeric_required.add_argument("-fe",type=str, help="Specify the Energy file", metavar="FILE", required=True)
numeric_required.add_argument("-fr",type=str, help="Specify the NACT Rho file", metavar="FILE", required=True)
numeric_required.add_argument("-fp",type=str, help="Specify the NACT Phi file", metavar="FILE", required=True)
numeric.add_argument("-h5", action='store_true', help="Write results in a HDF5 file (fast and efficient). (default behaviour)" )
numeric.add_argument("-txt" ,action="store_true", help="Write results in a text file.")
numeric.set_defaults(h5=False,txt=False)




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


#generation of analytic expression according to command line arguments
if args.choice =="ana":
    state = args.n 
    path  = args.p
    logger = make_logger("ADT Analytical program")
    try:
        adt_analytical(state, path, logger)
        print("Log saved in 'ADT.log'.")
    except Exception as e:
        logger.error("Program failed\n"+"-"*121)
        print("Program failed. %s"%e)



#calculation of ADT quantities according to command line arguments
if args.choice == "num":
    path = args.p 
    enrf = args.fe 
    rhof = args.fr 
    phif = args.fp
    outfile = args.o.strip("'")
    h5 = args.h5
    txt = args.txt
    if (h5==False and txt== False) : h5=True 

    logger = make_logger("ADT Numerical Program")
    try:
        adt_numerical(enrf, rhof, phif, path, outfile, logger, h5, txt)
        print("Log saved in 'ADT.log'.")
    except Exception as e:
        logger.error("Program failed\n"+"-"*121)
        print("Program failed. %s"%e)

