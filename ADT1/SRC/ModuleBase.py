#!/bin/python

######################################################################################################################################################
#This python program is written by Soumya Mukherjee, Bijit Mukherjee, Saikat Mukherjee, Subhankar Sardar and Satrajit Adhikari (corresponding author)
######################################################################################################################################################

################################################################################################################
def adiabatic(N):

    ''' This definition is built to generate a matrix form of adiabatic potential energy. It takes number of 
        electronic state (N) as an argument and returns the adiabatic potential energy (diagonal) matrix.'''

    mat = []

    for i in range(N):
      row = ()
      for j in range(N):
        if i == j:
          a = 'U(%r)' % (i+1)
          row += (a,)
        else:
          row += ('0',)

      mat += [(row)]
 
    return mat
   
 
################################################################################################################
def diabatic(mat):

    '''In this definition, the adiabatic potential energy matrix is similarity transformed to the diabatic 
       potential energy matrix. It accepts ADT matrix as input and returns the diabatic potential energy
       matrix.'''

    N = len(mat)

    umat = adiabatic(N)
    atmat = transpose(mat)

    w1 = multiply(atmat,umat)
    diabatic = multiply(w1,mat)

    return diabatic

 
################################################################################################################
def diff(elem,N):

    '''During implementing ADT condition, elements of ADT matrix are differentiated to evaluate the expression 
       of gradient of angles. This definition can be used to perform the required differentiation for N 
       dimensional sub-Hilbert space. This routine returns the differentiated form of any element of 
       ADT matrix.'''

    lltot = elem.split('+')

    convertlist=indextostate(N)   

    element = []

    for i in lltot:
      llsub = i.split('*')
      difflist = []

      for j in llsub:
        ll = j.split('(')
        ll1 = ll[-1].split(')')
        index = convertlist[int(ll1[0])]
        difflist.append(diffelem(j,index))

      totelem = []
      for k in range(len(difflist)):
        gradelem = []
        if k == 0:
          gradelem.append(difflist[k])
          for k1 in range(1,len(llsub)):
            gradelem.append(llsub[k1])
        else:
          for k1 in range(k):
            if k1 == 0:
              gradelem.append(llsub[k1])
            else:
              gradelem.append(llsub[k1])
          gradelem.append(difflist[k])
          for k2 in range(k+1,len(llsub)):
            gradelem.append(llsub[k2])        

        totelem.append('*'.join(gradelem))
           
      element.append('+'.join(totelem))
  
    elemfinal = '+'.join(element) 

    return elemfinal


################################################################################################################
def diffelem(trig,index):

    '''It takes a trigonometric function as argument and returns its derivative.'''    

    if 's' in trig:
      difftrig = trig.replace('s','c')
      difftrig += '*' + 'G(%s)' % index
    elif '-c' in trig:
      difftrig = trig.replace('-c','s')
      difftrig += '*' + 'G(%s)' % index
    elif 'c' in trig:
      difftrig = trig.replace('c','-s')  
      difftrig += '*' + 'G(%s)' % index

    return difftrig


################################################################################################################
def elemgradselect(mat):

    '''Although implementation of ADT condition results into a differential matrix equation, we can obtain only 
       ^{N}C_{2} number of independent equations (not N^{2}) by comparing the matrix elements. This routine 
       gathers necessary matrix elements, takes their derivative using 'diff' definition and finally collects 
       the outcome of the differentiation in a dictionary. That dictionary is returned at the end of such 
       processes.'''

    ll = {}
    N = len(mat)

    for i in range(N-1):
      a = 'G(%r,%r)' % (i+1,N)
      ll[a] = diff(mat[N-1][i],N)

    for i in range(N-2):
      for j in range(N-2,i,-1):
        a = 'G(%r,%r)' % (i+1,j+1)
        ll[a] = diff(mat[j][i],N)
 
    return ll


################################################################################################################
def elemselect(mat):
 
    '''This definition returns the list of lower triangle elements of any matrix (mat_{ij}, where i >= j).'''
 
    ll = []
    N = len(mat)

    for i in range(1,N):
      for j in range(i):
        ll.append(mat[i][j])

    return ll


################################################################################################################
def elemsum(elem1,elem2):

    '''This definition is designed to return the product of any two expressions.'''
    
    l1 = elem1.split('+')
    l2 = elem2.split('+')

    elem = []

    if elem1 == '0' or elem2 == '0':
      elemfinal = '0'
    else:
      for i in l1:
        for j in l2:
          if i == '1':
            elem.append(j)
          elif j == '1':
            elem.append(i)
          else:
            elem.append(i + '*' + j)
      elemfinal = '+'.join(elem) 

    return elemfinal
           

################################################################################################################
def elemtauselect(mat):

    '''Although implementation of ADT condition results into a differential matrix equation, we can obtain 
       only ^{N}C_{2} number of independent equations (not N^{2}) by comparing the matrix elements. This 
       definition fetches necessary matrix elements (product matrix of negative NACM and ADT matrix) and 
       store them in a dictionary. At the end of these processes, that dictionary is returned.'''

    ll = {}
    N = len(mat)
 
    for i in range(N-1):
      a = 'G(%r,%r)' % (i+1,N)
      ll[a] = mat[N-1][i]

    for i in range(N-2):
      for j in range(N-2,i,-1):
        a = 'G(%r,%r)' % (i+1,j+1)
        ll[a] = mat[j][i]
 
    return ll


################################################################################################################
def equation_complete(lhs,rhs,N):

    '''This definition is built to return the fully substituted forms of the ADT equations, but the major 
       demerit of this definition is high memory requirements, while constructing ADT equations for higher 
       dimensional sub-Hilbert spaces.'''   

    fl = open('ADT_EQUATIONS_COMPLETE.DAT','a')

    ll = []
    lldict = {}
    for i in range(1,N):
      for j in range(N,i,-1):
        leqn = []
        a = 'G(%r,%r)' % (i,j)
        multi = inver(lhs[a],a)
        rhsmod = elemsum(multi,rhs[a])
        leqn.append(rhsmod)
        lhselem = lhs[a]
        llnew = lhselem.split('+')
        for k in llnew:
          if a not in k:
            left = '(-1)' + '*' + multi + '*' + k
            leqn.append(left)
        sol = '+'.join(leqn)
        if ll:
          for b in ll:
            eqn = '(' + lldict[b] + ')'
            count = sol.count(b)
            for cindex in range(count):
              sol = sol.replace(b,eqn)
        fl.write('\n')
        fl.write(' GRAD_%r_%r = \n' % (i,j))
        length = len(sol)
        line = length/140
        for k1 in range(1,line+1):
          start = 0 + (k1-1)*140
          end = k1*140 
          segment = sol[start:end]
          fl.write(' ')
          fl.write(segment)
          fl.write('\n')
        fl.write(' ')
        fl.write(sol[line*140:])
        fl.write('\n')
        fl.write('\n')
   
        ll.append(a)
        lldict[a] = sol

    fl.close()


################################################################################################################
def equation_partial(lhs,rhs,N):

    '''This definition also returns the ADT equations in partially substituted forms, where the gradients in 
       the right hand side are not substituted by their expressions.'''   

    fl = open('ADT_EQUATIONS_PARTIAL.DAT','a')

    ll = []
    for i in range(1,N):
      for j in range(N,i,-1):
        leqn = []
        a = 'G(%r,%r)' % (i,j)
        multi = inver(lhs[a],a)
        rhsmod = elemsum(multi,rhs[a])
        leqn.append(rhsmod)
        lhselem = lhs[a]
        if ll:
          for b in ll:
            eqn = 'Sol(%r)' % b
            count = lhselem.count(b)
            for cindex in range(count):
              lhselem = lhselem.replace(b,eqn)
          llnew = lhselem.split('+')
          for k in llnew:
            if a not in k:
              left = '(-1)' + '*' + multi + '*' + k
              leqn.append(left)
        sol = '+'.join(leqn)
        fl.write('\n')
        fl.write(' GRAD_%r_%r = \n' % (i,j))
        length = len(sol)
        line = length/140
        for k1 in range(1,line+1):
          start = 0 + (k1-1)*140
          end = k1*140 
          segment = sol[start:end]
          fl.write(' ')
          fl.write(segment)
          fl.write('\n')
        fl.write(' ')
        fl.write(sol[line*140:])
        fl.write('\n')
        fl.write('\n')
   
        ll.append(a)

    fl.close()


################################################################################################################
def extractgrad(eqn,N,gradrow):

    '''This definition gives the coefficient matrix of the gradient of angles (ADT angles).'''

    fl = open('GRADCOEFF.DAT','a')

    ll = eqn.split('+')
    ND = N*(N-1)/2
    gradcolumn = 0

    for i in range(2,N+1):
      for j in range(1,i):
        gradcolumn += 1
        a = 'G(%r,%r)' % (j,i)
        cinter = []
        for k in range(len(ll)):
          lsplit = ll[k].split('*')
          if a in lsplit:
            lsplit.remove(a)
            E = '*'.join(lsplit)
            cinter.append(E)
        if cinter:
          coeff = '+'.join(cinter)
        else:
          coeff = '0'
        fl.write('\n')
        fl.write(' GC_%r_%r = \n' % (gradrow,gradcolumn))
        length = len(coeff)
        line = length/140
        for k1 in range(1,line+1):
          start = 0 + (k1-1)*140
          end = k1*140 
          segment = coeff[start:end]
          fl.write(' ')
          fl.write(segment)
          fl.write('\n')
        fl.write(' ')
        fl.write(coeff[line*140:])
        fl.write('\n')
        fl.write('\n')
    
    fl.close()


################################################################################################################
def extracttau(eqn,N,taurow):

    '''This definition gives the coefficient matrix of the nonadiabatic coupling terms (NACTs).'''

    fl = open('TAUCOEFF.DAT','a')

    ll = eqn.split('+')
    ND = N*(N-1)/2
    taucolumn = 0

    for i in range(2,N+1):
      for j in range(1,i):
        taucolumn += 1
        a = '(TAU(%r,%r))' % (j,i)
        cinter = []
        for k in range(len(ll)):
          lsplit = ll[k].split('*')
          if a in lsplit:
            lsplit.remove(a)
            E = '*'.join(lsplit)
            cinter.append(E)
        if cinter:
          coeff = '+'.join(cinter)
        else:
          coeff = '0'
        fl.write('\n')
        fl.write(' TC_%r_%r = \n' % (taurow,taucolumn))
        length = len(coeff)
        line = length/140
        for k1 in range(1,line+1):
          start = 0 + (k1-1)*140
          end = k1*140 
          segment = coeff[start:end]
          fl.write(' ')
          fl.write(segment)
          fl.write('\n')
        fl.write(' ')
        fl.write(coeff[line*140:])
        fl.write('\n')
        fl.write('\n')
    
    fl.close()


################################################################################################################
def indextostate(N):

    '''Overall ADT matrix for any arbitrary number of coupled electronic states is generated by multiplying 
       the elementary rotation matrices in a definite order. Their position in that multiplication order is 
       clearly dictated by the electronic states involved in the corresponding mixing angle. As an example, 
       for three electronic states, the order of multiplication is A_{12}.A_{13}.A_{23}. In order to 
       correlate the order of appearance of the elementary ADT matrix and the involved electronic states, 
       a dictionary is returned by this definition. [For A_{12}.A_{13}.A_{23}, the returned dictionary 
       will be {1:1,2,  2:1,3,  3:2,3}]'''

    indict = {}

    counter = 0
    for i in range(2,N+1):
      for j in range(1,i):
        a = '%r,%r' % (j,i)
        counter += 1
        indict[counter]=a

    return indict

                
################################################################################################################
def inver(element,grad):

    '''This definition acquires an expression (grad), replaces a trigonometric function (element) with its 
       inverse and finally, returns the modified version of 'grad'.'''

    ll = element.split('+')
    lnew = []
    for i in ll:
      if grad in i:
        llsub = i.split('*')
        for e in llsub:
          if 'c' in e:
            inew = e.replace('c','ic')
            lnew.append(inew)
          elif 's' in e:
            inew = e.replace('s','is')
            lnew.append(inew)

    eleminver = '*'.join(lnew)

    return eleminver

                                                                                             
################################################################################################################
def matelem(mat1,mat2,i,j):

    '''This definition returns (i,j) element of product of two matrices, namely, 'mat1' and 'mat2'.'''
 
    N = len(mat1)
    elem = []

    for k in range(N):
      elem.append(elemsum(mat1[i][k],mat2[k][j]))
 
    for n in range(elem.count('0')):
      elem.remove('0') 

    if not elem:
      elemnew = '0'
    else:
      elemnew = '+'.join(elem)

    return elemnew  


################################################################################################################
def matman(N):

    '''This produces the multiplied form of ADT matrix for any arbitrary number of coupled electronic 
       states.'''

    counter = 0
    mat1 = unitmat(N)
    for i in range(2,N+1):
      for j in range(1,i):
        counter += 1
        mat2 = matrix(j,i,N,counter) 
        lmat = multiply(mat1,mat2)
        mat1 = lmat
   
    return lmat

  
################################################################################################################
def matrix(first,second,N,counter):

    '''This definition generates elementary rotation matrices, which are multiplied to construct the overall
       ADT matrix.'''
     
    mat = []
    for i in range(1,N+1):
      row = ()
      for j in range(1,N+1):
        if (i == first and j == first) or (i == second and j == second):
          a = '(c(%r))' % counter
          row += (a,)  
        elif i == first and j == second:
          a = '(s(%r))' % counter
          row += (a,) 
        elif i == second and j == first:
          a = '(-s(%r))' % counter
          row += (a,)
        elif i == j:
          row += ('1',)
        else:
          row += ('0',)
      mat += [(row)]
  
    return mat
 
 
################################################################################################################
def multiply(mat1,mat2):

    '''This definition returns multiplied expression of two matrices ('mat1' and 'mat2').'''

    N = len(mat1)
    newmat = []

    for i in range(N):
      matrow = ()
      for j in range(N):
        matrow += (matelem(mat1,mat2,i,j),)
      
      newmat += [(matrow)]        
   
    return newmat


################################################################################################################
def nacm(N):

   '''This returns nonadiabatic coupling matrix (NACM) for N electronic state sub-Hilbert space.'''

   ll = []
   
   for i in range(N):
     row = ()
     for j in range(N):
       if i == j:
         row += ('0',)
       elif i > j:
         a = '(-1)*(TAU(%r,%r))' % (j+1,i+1)
         row += (a,)
       else:
         a = '(TAU(%r,%r))' % (i+1,j+1)
         row += (a,)

     ll += [(row)]

   return ll

  
################################################################################################################
def negative(mat):

    '''This definition is implemented to generate a new matrix, where the input matrix is multiplied with 
       (-1).'''

    N = len(mat)
    matnew = []
    
    for i in range(N):
      row = ()
      for j in range(N):
        a = elemsum('(-1)',mat[i][j])
        row += (a,)
 
      matnew += [(row)]

    return matnew    
  

################################################################################################################
def transpose(mat):

    '''This definition returns the transpose of a matrix (mat).'''

    matconj = []

    for i in range(len(mat)):
      row = ()
      for j in range(len(mat)):
        if i == j:
          row += (mat[i][j],)
        else:
          row += (mat[j][i],)

      matconj += [(row)]

    return matconj


################################################################################################################
def unitmat(N):

    '''Unit matrix of dimension N can be obtained by this definition.'''

    matfirst = []
    for i in range(N):
      row = ()
      for j in range(N):
        if i == j:
          row += ('1',)
        else:
          row += ('0',)
 
      matfirst += [(row)]

    return matfirst


################################################################################################################


