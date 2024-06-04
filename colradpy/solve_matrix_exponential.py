import numpy as np

def solve_matrix_exponential(matrix,td_n0,td_t):
    """This definition will solve a 4 dimensional matrix using the matrix expoentiaiation method
       

       R. LeVeque, Finite Difference Methods for Ordinary and Par-
       tial Differential Equations: Steady-State and Time-Dependent
       Problems (Classics in Applied Mathematics Classics in Applied
       Mathemat), Society for Industrial and Applied Mathematics,
       Philadelphia, PA, USA, 2007.



    Args:
      :param matrix: The 4d matrix to be solved
      :type matrix: 4d matrix x,y,temp,dens

      :param td_n0: The initial fractional populations
      :type td_n0: float array

      :param td_t: array of times for the solution
      :type metas: float array


    Returns:
      This returns three arrays the time dependent populations, eigenvals and eigenvectorsimport matplotlib.pyplot as plt

    """
    #calculate the eigenvectors and values for the matrix
    #the last two dimensions have to be the ij parts so we need to transpose matrix
    #string logic for the einsums based on matrix dimensions        
    if(len(np.shape(matrix)) ==4):
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,3,0,1))

        eval_str = 'klj'
        evec_str = 'klij'
        vt_str1 = 'klj'
        pop_str='ikl'
    if(len(np.shape(matrix)) ==3):
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,0,1))

        eval_str = 'lj'
        evec_str = 'lij'
        vt_str1 = 'lj'
        pop_str = 'il'

    if(len(np.shape(td_t)) > 1):    
        if(np.shape(td_t)[0] ==1):
            td_t = td_t[0,:]
            td_t_str = 'l'
            vt_str2 = ''
        else:
            td_t_str = 'lt'
            pop_str = pop_str[0:1] +'t' + pop_str[1:]
            vt_str2='t'
    else:
        td_t_str = 't'
        vt_str2 = 't'
        pop_str = pop_str[0:1] +'t' + pop_str[1:]

    v0 = np.dot(np.linalg.inv(eigenvectors),td_n0)        
    #vt = v0[:,:,:,None]*np.exp(eigenvals[:,:,:,None]*td_t)
    vt = np.einsum(eval_str+ ',' +vt_str1+vt_str2 +'->'+vt_str1+vt_str2,v0,
                   np.exp(np.einsum(eval_str+','+td_t_str+'->'+vt_str1+vt_str2  ,eigenvals,td_t)))

    td_pop = np.einsum(evec_str+','+vt_str1+vt_str2+'->'+pop_str  , eigenvectors,vt)
    
    return td_pop, eigenvals, eigenvectors


def eval_matrix_exponential_solution(time, n0, eigenvalues, eigenvectors):
    v0 = np.dot(np.linalg.inv(eigenvectors), n0)
    vt = np.einsum('klj,kljt->kljt', v0, np.exp(np.einsum('klj,t->kljt', eigenvalues, time)))
    pops = np.einsum('klij,kljt->itkl', eigenvectors, vt)
    pops[pops < 0] = 0
    return pops


def solve_matrix_exponential_steady_state(matrix):
    """ This definition will solve the steady state solutions of a 4d matrix
       

       R. LeVeque, Finite Difference Methods for Ordinary and Par-
       tial Differential Equations: Steady-State and Time-Dependent
       Problems (Classics in Applied Mathematics Classics in Applied
       Mathemat), Society for Industrial and Applied Mathematics,
       Philadelphia, PA, USA, 2007.

    In steady state only the smallest eigenval survives
    so find the smallest eigenval and the corresponding eigenvector 
    then solve for the populations with just that eigen value.



    Args:
      :param matrix: The 4d matrix to be solved
      :type matrix: 4d matrix[x,y,temp,dens]

    Returns:
      This returns three arrays the steady state populations, eigenvals and eigenvectors

    """

    #becuase this is steady state the initial populations don't matter
    td_n0 = np.zeros(np.shape(matrix)[0]) 
    td_n0[0] = 1.
    #sort on the eigenval index to find the longest lived state
    #that will give the equilibrium populuations
    if(len(np.shape(matrix)) ==4):
        
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,3,0,1))
        axis=2
    if(len(np.shape(matrix)) ==3):
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,0,1))
        axis=1

        
    v0 = np.dot(np.linalg.inv(eigenvectors),td_n0)

    index = list(np.ix_(*[np.arange(i) for i in eigenvals.shape]))
    index[axis] = np.abs(eigenvals).argsort(axis)
    index = tuple(index)  # Needs to be tuple for multidimensional indexing to work

    if(len(np.shape(matrix)) ==4):    
        ev = eigenvectors.transpose(0,1,3,2)[index]#egienvectors sorted on eigenvals
        ss_pop = np.einsum('kl,klj->klj',v0[index][:,:,0],ev[:,:,0,:]).transpose(2,0,1)
    if(len(np.shape(matrix)) ==3):
        ev = eigenvectors.transpose(0,2,1)[index]#egienvectors sorted on eigenvals
        ss_pop = np.einsum('k,kj->kj',v0[index][:,0],ev[:,0,:]).transpose(1,0)

    return ss_pop, eigenvals, eigenvectors






def solve_matrix_exponential_source(matrix, td_n0, source, td_t):
    """This definition will solve a 4 dimensional matrix using the matrix expoentiaiation method
       when a source term is also included
       
       This is a slight modification to R. LeVeque 2007, see Johnson thesis


    Args:
      :param matrix: The 4d matrix to be solved
      :type matrix: 4d matrix x,y,temp,dens

      :param td_n0: The initial fractional populations
      :type td_n0: float array

      :param source: The source of particles into the different states.
      :type source: float array


      :param td_t: array of times for the solution
      :type metas: float array


    Returns:
      This returns three arrays the time dependent populations, eigenvals and eigenvectors

    """
    
    #calculate the eigenvectors and values for the matrix
    #the last two dimensions have to be the ij parts so we need to transpose matrix
    #string logic for the einsums based on matrix dimensions        
    if(len(np.shape(matrix)) ==4):
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,3,0,1))

        eval_str = 'klj'
        evec_str = 'klij'
        vt_str1 = 'klj'
        pop_str='ikl'
    if(len(np.shape(matrix)) ==3):
        eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,0,1))

        eval_str = 'lj'
        evec_str = 'lij'
        vt_str1 = 'lj'
        pop_str = 'il'

    if(len(np.shape(td_t)) > 1):    
        if(np.shape(td_t)[0] ==1):
            td_t = td_t[0,:]
            td_t_str = 'l'
            vt_str2 = ''
        else:
            td_t_str = 'lt'
            pop_str = pop_str[0:1] +'t' + pop_str[1:]
            vt_str2='t'
    else:
        td_t_str = 't'
        vt_str2 = 't'
        pop_str = pop_str[0:1] +'t' + pop_str[1:]

    #the math equations for this solution
    #are in Johnson thesis 2020 apprendix A
    CC = np.dot(np.linalg.inv(eigenvectors),source)
    V0 = np.dot(np.linalg.inv(eigenvectors),td_n0)

    v_tmp = []
    
    #this number might need to change to reject larger eigen values
    #there are numerical problems with small eigenvalues so just treat them like zero eigenvalues
    eig_zero_ind = np.where(np.abs(eigenvals) < 1.e-10)#np.finfo(np.float64).eps*40000)
    eig_non_zero = np.copy(eigenvals)
    eig_non_zero[eig_zero_ind] = 1.

    amplitude_zer = np.zeros_like(V0)
    amplitude_zer[eig_zero_ind] =V0[eig_zero_ind]
    
    CC_zer = np.zeros_like(CC)
    CC_zer[eig_zero_ind] = CC[eig_zero_ind]
    V0[eig_zero_ind] = 0.
    CC[eig_zero_ind] = 0.
    
    amplitude_non = V0 +CC/eig_non_zero

    ttmp = CC/eig_non_zero
    
    #these are the coded up equations but allowing for the different combinations of matrix and
    # time dimensions
    # v_non = amplitude_non[:,:,:,None]*np.exp(eig_non_zero[:,:,:,None]*td_t) - \
    #                           np.delete(CC,eig_zero_ind,axis=2)[:,:,:,None]/eig_non_zero[:,:,:,None]
    #v_zer = CC[:,:,eig_zero_ind[2]][:,:,:,None]*td_t + amplitude_zer[:,:,:,None]

    #if there is a time dependence the ttmp array needs to add on another dimension
    if(vt_str2 == 't'):

        v_non = np.einsum( eval_str+ ',' +vt_str1+vt_str2 +'->'+vt_str1+vt_str2, amplitude_non,
                       np.exp(np.einsum(eval_str+','+td_t_str+'->'+vt_str1+vt_str2  ,eig_non_zero,td_t))) - \
                ttmp[...,None]
    else:
        v_non = np.einsum( eval_str+ ',' +vt_str1+vt_str2 +'->'+vt_str1+vt_str2, amplitude_non,
                       np.exp(np.einsum(eval_str+','+td_t_str+'->'+vt_str1+vt_str2  ,eig_non_zero,td_t))) - \
                ttmp

    #if there is a time dependence the amplitude_zer array needs to add on another dimension    
    if(eig_zero_ind[0].size > 0):

        if(vt_str2 == 't'):        
            v_zer = np.einsum(eval_str+','+td_t_str+'->'+vt_str1+vt_str2, CC_zer,td_t) + amplitude_zer[...,None]
        else:
            v_zer = np.einsum(eval_str+','+td_t_str+'->'+vt_str1+vt_str2, CC_zer,td_t) + amplitude_zer
            
        v_non[eig_zero_ind] = v_zer[eig_zero_ind]

    v_non[np.abs(v_non) <1.e-10] = 0

    td_pop = np.einsum(evec_str+','+vt_str1+vt_str2+'->'+pop_str  , eigenvectors,v_non)
    
    return td_pop, eigenvals,eigenvectors
