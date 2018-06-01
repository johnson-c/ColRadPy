import numpy as np
from scipy.linalg import qr
import scipy
from scipy import linalg, matrix

def qr_null(A, tol=None):
    Q, R, P = qr(A.T, mode='full', pivoting=True)
    tol = np.finfo(R.dtype).eps if tol is None else tol
    rnk = min(A.shape) - np.abs(np.diag(R))[::-1].searchsorted(tol)
    return Q[:, rnk:].conj()
def null(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)
def null(A, eps=1e-1):
    u, s, vh = scipy.linalg.svd(A)
    padding = max(0,np.shape(A)[1]-np.shape(s)[0])
    null_mask = np.concatenate(((s <= eps), np.ones((padding,),dtype=bool)),axis=0)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)


def null(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns





'''
charge_state_files = np.array(['/home/curtis/spirit/adf04-10apr17_mod'])
charge_state_files = np.array(['/home/curtis/spirit/wdarc7_mod'])
charge_state_files = np.array(['/home/curtis/spirit/wlike_mons11#w0.dat'])
charge_state_files =  '/home/curtis/Downloads/adf04_WI_16.10.17_ecip'
charge_state_files = '/home/curtis/spirit/wdx/tdcc/adf04_WI_5.12-3.17_tdcc'

charge_state_files = '/home/curtis/spirit/mo_i_figure_390nm_connors_paper/adf04'
#np.array(['/home/curtis/lewis/he_0/mom97_ls#he0.dat','/home/curtis/lewis/he_0/mom97_ls#he1.dat']) 


#sig = np.array([0,1])
#sig = np.linspace(0,39,40,dtype='int')
sig = np.array([0,1])
nsigmas = np.array([len(sig)])
gcrs = []

temperature_grid = np.array([6])#temp_arr
electron_den = np.array([1.e11]) #dens_arr

#for i in range(0,len(charge_state_files)):

'''

def ioniz(gc):
    gcrs=[]
    nsigmas = np.array([len(gc['metas'])])
    sig = gc['metas']
    gcrs.append(gc)
    temperature_grid = gc['user_temp_grid']
    electron_den = gc['user_dens_grid']

    ionbal = np.zeros((np.sum(nsigmas)+1,np.sum(nsigmas)+1,
                       len(gcrs[0]['user_temp_grid']),
                       len(gcrs[0]['user_dens_grid'])))

    '''
    i = 0
    for k in range(0, np.sum(nsigmas)):
        if(k > np.sum(nsigmas[0:i])):
            i = i + 1
        for j in range(0,np.sum(nsigmas)):
            print(k,j,i)

            if(k != j):
                if(gcrs[i]['qcd'].any() and j <len(gcrs[i]['qcd']) ):
                    print('here')
                    ionbal[k,j,:,:] = gcrs[i]['qcd'][j,k,:,:]
            else:

            ionbal[k,j,:,:] = -1*np.sum(gcrs[i]['scd'][j,:,:,:],axis=1)
            if(gcrs[i]['qcd'].any()):
                ionbal[k,j,:,:] = ionbal[k,j,:,:] + -1*np.sum(gcrs[i]['qcd'][j,:,:,:],axis=1)
    '''
    i = 0
    p = 0
    q = 0
    r = 0
    s = 0
    t = 0

    for k in range(0, np.sum(nsigmas)):
        if(k > np.sum(nsigmas[0:i+1]-1)):
            i = i + 1
            p = 0
            q = 0
            r = 0
        for j in range(0,np.sum(nsigmas)):

            #put the -qcd on the diagonal
            if( k==j and j < np.sum(nsigmas) and gcrs[i]['qcd'].any()):
                ionbal[k,j,:,:] = ionbal[k,j,:,:] + -1*(np.sum(gcrs[i]['qcd'][p,:,:,:],axis=0))
                p = p + 1

            if( k!=j and gcrs[i]['qcd'].any()):
                if(k<nsigmas[i] and j <nsigmas[i]):
                    if(k>j):
                        ionbal[k,j,:,:] = ionbal[k,j,:,:] + gcrs[i]['qcd'][j,k-1,:,:]
                    else:
                        ionbal[k,j,:,:] = ionbal[k,j,:,:] + gcrs[i]['qcd'][j,k,:,:]                    

            if( k!=j and gcrs[s]['scd'].any() and (k>np.sum(nsigmas[0:s+1])-1 or j >np.sum(nsigmas[0:s+1])-1) ):
                if(k < j):
                    ionbal[k,j,:,:] = ionbal[k,j,:,:] + gcrs[s]['acd'][k,j-nsigmas[s],:,:]
                if(k > j):
                    ionbal[k,j,:,:] = ionbal[k,j,:,:] + gcrs[s]['scd'][j,k-nsigmas[s],:,:]

            #put the -scd on the diagonal
            if( k==j and gcrs[i]['scd'].any()):
                ionbal[k,j,:,:] = ionbal[k,j,:,:] + -1*(np.sum(gcrs[i]['scd'][q,:,:,:],axis=0))
                q = q + 1

            #put the -acd on the diagonal
            if( k==j and gcrs[i]['acd'].any() and k >= nsigmas[0]-1):

                ionbal[np.sum(nsigmas[0:i+1]),np.sum(nsigmas[0:i+1]),:,:] = ionbal[np.sum(nsigmas[0:i+1]),np.sum(nsigmas[0:i+1]),:,:]\
                -1*(np.sum(np.sum(gcrs[i]['acd'][:,:,:,:],axis=0),axis=0))

                r = r+ 1

    u = 0            
    for k in range( np.sum(nsigmas) - nsigmas[-1]  ,np.sum(nsigmas)):
        ionbal[-1,k,:,:] = ionbal[-1,k,:,:] + gcrs[-1]['scd'][u,0,:,:]
        ionbal[k,-1,:,:] = ionbal[k,-1,:,:] + gcrs[-1]['acd'][u,0,:,:]

        u = u+1



    #include source term


    

    t_min = 0
    t_max = .001
    dt = (np.max(ionbal))
    t_steps = np.array([0,.001,t_max])
    t_steps = np.linspace(0,.0001,1000)
    td_pop = np.zeros( len(ionbal[0]))
    td_pop[0] = 1.
    td_pop_correct = np.zeros( (len(td_pop),len(t_steps),len(temperature_grid),len(electron_den)))
    for e in range(0,len(electron_den)):
        for t in range(0,len(temperature_grid)):
            eigenval,eigenvectors = np.linalg.eig(ionbal[:,:,t,e]*electron_den)
            v0 = np.dot(np.linalg.inv(eigenvectors),td_pop)
            vt = v0[:,None]*np.exp(eigenval[:,None]*t_steps)
            td_pop_correct[:,:,t,e] = np.dot(eigenvectors,vt)
    return t_steps,td_pop_correct
    
    ##########################################################################################
    #
    # time depepent part
    #
    ##########################################################################################

    t_min = 0
    t_max = .001
    dt = (np.max(ionbal))*100.
    print(dt)
    nsteps = (t_max - t_min) /dt
    nsteps = nsteps
    t_steps = np.linspace(t_min,t_max,nsteps)

    td_pop = np.zeros( (np.shape(t_steps) + np.shape(ionbal[0]) ))
    td_pop[0,0,:,:] = 1

    for t in range(0,len(gcrs[0]['user_temp_grid'])):
        for n in range(0,len(gcrs[0]['user_dens_grid'])):
            for i in range(1, len(t_steps)):
                td_pop[i,:,t,n] = td_pop[i-1,:,t,n] + (1/dt)*np.dot(ionbal[:,:,t,n],td_pop[i-1,:,t,n])








                


    ##########################################################################################
    #
    # time depepent part with source
    #
    ##########################################################################################

    source =1.

    td_pop_s = np.zeros( (np.shape(t_steps) + np.shape(ionbal[0]) ))

    source_arr = np.zeros(len(sig)+1)
    source_arr[0] = source
    td_pop_s[0,0,:,:] = 1.
    for t in range(0,len(gcrs[0]['user_temp_grid'])):
        for n in range(0,len(gcrs[0]['user_dens_grid'])):
            for i in range(1, len(t_steps)):
                td_pop_s[i,:,t,n] = td_pop_s[i-1,:,t,n] + (1/dt)*np.dot(ionbal[:,:,t,n],td_pop_s[i-1,:,t,n]) + source_arr 
    return t_steps,td_pop,td_pop_s,td_pop_correct
