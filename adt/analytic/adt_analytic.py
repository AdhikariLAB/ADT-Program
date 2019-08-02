from __future__ import absolute_import, unicode_literals, division, print_function
__doc__='''

This python script is specifically implemented for devising analytic expressions of eight adiabatic to diabatic
transformation (ADT) quantities for any 'N' number of coupled electronic states.

Short description of eight definitions:

adt1 :  gives elements of adiabatic potential energy matrix
adt2 :  gives elements of nonadiabatic coupling matrix (NACM)
adt3 :  gives elements of adiabatic to diabatic transformation (ADT) matrix
adt4 :  gives partially substituted forms of ADT equations
adt5 :  gives complete forms of ADT equations
adt6 :  gives elements of coefficient matrix of gradient of ADT angles
adt7 :  gives elements of coefficient matrix of nonadiabatic coupling terms (NACTs)
adt8 :  gives elements of diabatic potential energy matrix

'''
__authors__  = '''
Koushik Naskar, Soumya Mukherjee, Bijit Mukherjee, Satyam Ravi, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari
'''


from adt.analytic import anamod
from time import time


#Main switcher function

# def adt_analytical(N,p, logger):
#     start = time()

#     {1: adt1, 2: adt2, 3: adt3, 4: adt4,
#      5: adt5, 6: adt6, 7: adt7, 8: adt8
#     }[p](N,logger)

#     logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)


# 6 available options, removing first two from earlier
def adt_analytical(N,p, logger):
    start = time()

    {1: adt3, 2: adt4,
     3: adt5, 4: adt6, 5: adt7, 6: adt8
    }[p](N,logger)

    logger.info("Program completed successfully in %.5f seconds\n"%(time()-start)+"-"*121)


#This definition returns elements of adiabatic potential energy matrix

# def adt1(N,logger):
#     logger.info("Deriving adiabatic potential energy matrix")
#     U_Matrix = anamod.adiabatic(N)

#     txt =""
#     for i in range(1,N+1):
#       for j in range(1,N+1):
#         txt += ' U_%s_%s = ' % (i,j)
#         txt += U_Matrix[i-1][j-1] +"\n\n"

#     logger.info("Writing results in 'U_MATRIX.DAT'")

#     with open('U_MATRIX.DAT', "w") as f:
#         f.write(txt)


# #This definition returns elements of nonadiabatic coupling matrix (NACM)

# def adt2(N,logger):
#     logger.info("Deriving NACM")
#     NACM = anamod.nacm(N)

#     txt =""
#     for i in range(1,N+1):
#       for j in range(1,N+1):
#         txt += ' NACM_%s_%s = ' % (i,j)
#         txt += NACM[i-1][j-1] + "\n\n"

#     logger.info("Writing results in 'NACM.DAT'")
#     with open('NACM.DAT','w') as f:
#         f.write(txt)


#This definition returns elements of adiabatic to diabatic transformation (ADT) matrix

def adt3(N,logger):
    logger.info("Deriving ADT matrix elements")
    txt   = ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
            .format(c=count, i1=index1, i2=index2)

    A_Matrix = anamod.matman(N)

    for i in range(1,N+1):
      for j in range(1,N+1):
        txt += ' A_%s_%s = \n' % (i,j)
        E = A_Matrix[i-1][j-1]
        txt += "\n".join(E[i:i+140] for i in range(0,len(E),140))+"\n\n"

    logger.info("Writing results in 'A_MATRIX.DAT'")
    with open('A_MATRIX.DAT','w') as f:
        f.write(txt)


#This definition returns partially substituted forms of ADT equations

def adt4(N,logger):
    logger.info("Deriving partially substituted ADT equations")
    A_Matrix = anamod.matman(N)

    #differentiating the relevant elements of ADT matrix
    LHSelems = anamod.elemgradselect(A_Matrix)

    #formation of nonadiabatic coupling matrix (NACM)
    NACM = anamod.nacm(N)

    #multiplication of negative of NACM and ADT matrix
    R_Matrix = anamod.multiply(anamod.negative(NACM),A_Matrix)

    #collecting the relevant elements of the above product matrix
    RHSelems = anamod.elemtauselect(R_Matrix)

    #writing the partially substituted form of ADT equations


    txt = ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt +="""
    c({c}) = cos(theta({i1},{i2}))

    s({c}) = sin(theta({i1},{i2}))

    ic({c})= sec(theta({i1},{i2}))

    is({c})= cosec(theta({i1},{i2}))
    """.format(c=count, i1=index1, i2=index2)


    with open("ADT_EQUATIONS_PARTIAL.DAT", "w") as f:
        f.write(txt)

    logger.info("Derivation in process... Writing results in 'ADT_EQUATIONS_PARTIAL.DAT'")
    Equations = anamod.equation_partial(LHSelems,RHSelems,N)


#This definition returns complete forms of ADT equations

def adt5(N,logger):
    logger.info("Deriving complete form of ADT equations")
    #formation of adiabatic to diabatic transformation (ADT) matrix
    A_Matrix = anamod.matman(N)

    #differentiating the relevant elements of ADT matrix
    LHSelems = anamod.elemgradselect(A_Matrix)

    #formation of nonadiabatic coupling matrix (NACM)
    NACM = anamod.nacm(N)

    #multiplication of negative of NACM and ADT matrix
    R_Matrix =anamod.multiply(anamod.negative(NACM),A_Matrix)

    #collecting the relevant elements of the above product matrix
    RHSelems = anamod.elemtauselect(R_Matrix)

    #writing the completely substituted form of ADT equations
    txt = ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt +="""
    c({c}) = cos(theta({i1},{i2}))

    s({c}) = sin(theta({i1},{i2}))

    ic({c})= sec(theta({i1},{i2}))

    is({c})= cosec(theta({i1},{i2}))
    """.format(c=count, i1=index1, i2=index2)

    with open("ADT_EQUATIONS_COMPLETE.DAT", "w") as f:
        f.write(txt)

    logger.info("Derivation in process... Writing results in 'ADT_EQUATIONS_COMPLETE.DAT'")
    Equations = anamod.equation_complete(LHSelems,RHSelems,N)


#This definition returns elements of coefficient matrix of gradient of ADT angles

def adt6(N,logger):
    logger.info("Deriving coefficient matrix of gradient of ADT angles")
    A_Matrix = anamod.matman(N)

    #collecting the relevant elements of ADT matrix
    Aelems = anamod.elemselect(A_Matrix)

    #differentiating the collected elements of ADT matrix

    LHSelems = [anamod.diff(el,N) for el in Aelems]

    #writing the coefficient matrix of the gradient of the ADT angles

    txt= ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
            .format(c=count, i1=index1, i2=index2)

    with open('GRADCOEFF.DAT','w') as f:
        f.write(txt)

    logger.info("Derivation in process... Writing results in 'GRADCOEFF.DAT'")
    for j,i in enumerate(LHSelems,1):
        anamod.extractgrad(i,N,j)


#This definition returns elements of coefficient matrix of nonadiabatic coupling terms (NACTs)

def adt7(N,logger):
    logger.info("Deriving elements of coefficient matrix of NACTs")
    A_Matrix = anamod.matman(N)

    #formation of nonadiabatic coupling matrix (NACM)
    NACM = anamod.nacm(N)

    #multiplication of negative of NACM and ADT matrix
    R_Matrix = anamod.multiply(anamod.negative(NACM),A_Matrix)

    #collecting the relevant elements of the above product matrix
    RHSelems = anamod.elemselect(R_Matrix)

    #writing the coefficient matrix of the nonadiabatic coupling terms (NACT)
    txt= ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
            .format(c=count, i1=index1, i2=index2)

    with open('TAUCOEFF.DAT','w') as f:
        f.write(txt)

    logger.info("Derivation in process... Writing results in 'TAUCOEFF.DAT'")
    for j,i in enumerate(RHSelems,1):
        anamod.extracttau(i,N,j)



#This definition returns elements of diabatic potential energy matrix

def adt8(N,logger):
    logger.info("Deriving elements of diabatic potential energy matrix")
    #formation of adiabatic potential energy matrix
    A_Matrix = anamod.matman(N)

    #formation of diabatic potential energy matrix
    W_Matrix = anamod.diabatic(A_Matrix)

    #writing of diabatic potential energy matrix elements
    fl = open('W_MATRIX.DAT','w')

    txt= ""
    count = 0
    for index1 in range(2,N+1):
      for index2 in range(1,index1):
        count += 1
        txt += 'c({c}) = cos(theta({i2},{i1}))\n\ns({c}) = sin(theta({i2},{i1}))\n\n'\
            .format(c=count, i1=index1, i2=index2)

    fl.write(txt)

    logger.info("Writing results in 'W_MATRIX.DAT'")
    for i in range(1,N+1):
        for j in range(1,N+1):
            txt = ' W_%s_%s = \n' % (i,j)
            E = W_Matrix[i-1][j-1]
            txt += "\n".join(E[i:i+140] for i in range(0,len(E),140)) + "\n\n"
            fl.write(txt)

########################################################################################################################
