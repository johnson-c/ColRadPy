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
      This returns three arrays the time dependent populations, eigenvals and eigenvectors

    """
    
    eigenvals, eigenvectors = np.linalg.eig(matrix.transpose(2,3,0,1))
    v0 = np.dot(np.linalg.inv(eigenvectors),td_n0)
    vt = v0[:,:,:,None]*np.exp(eigenvals[:,:,:,None]*td_t)
    td_pop = np.einsum('klij,kljt->itkl', eigenvectors, vt)
    eigenvals = eigenvals.transpose(2,0,1)
    eigenvectors = eigenvectors.transpose(2,3,0,1)
    return td_pop, eigenvals, eigenvectors


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
    
    eigenvals,eigenvectors = np.linalg.eig(matrix.transpose(2,3,0,1))

    CC = np.dot(np.linalg.inv(eigenvectors),source)
    V0 = np.dot(np.linalg.inv(eigenvectors),td_n0)

    eig_zero_ind = np.where(eigenvals == 0)            
    eig_non_zero = np.delete(eigenvals, eig_zero_ind, axis=2)

    amplitude_non = np.delete(V0,eig_zero_ind,axis=2) + np.delete(CC,eig_zero_ind,axis=2)/eig_non_zero
    amplitude_zer = V0[:,:,eig_zero_ind[2]]
    
    v_non = amplitude_non[:,:,:,None]*np.exp(eig_non_zero[:,:,:,None]*td_t) - \
                               np.delete(CC,eig_zero_ind,axis=2)[:,:,:,None]/eig_non_zero[:,:,:,None]
    v_zer = CC[:,:,eig_zero_ind[2]][:,:,:,None]*td_t + amplitude_zer[:,:,:,None]
    v = np.insert(v_non,eig_zero_ind[2],v_zer,axis=2)
    td_pop = np.einsum('klij,kljt->itkl', eigenvectors,v)
    return td_pop, eigenvals,eigenvectors
