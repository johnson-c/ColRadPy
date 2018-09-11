################################################################################
# file name         : colrad.py
# author            : Curt Johnson
# description       : This code solves the collsion radiative matrix
# version           : 0.3
# python version    : 2.7.12 ipython 2.4.1
# dependencies      : read_adf04.py, numpy, scipy interpolate, matplotlib
# license           : to be freed after code is stable and papers are published
#
# This code originally went through different file names and versions that
# incorprated different functionalitiy before it was under version control.
# These are included in the repos for historical and bug tracking reasons.
# The order went adf04_read.py, adf04_read1.py, adf04_dens_grid.py,
# adf04_dens_temp_grid_tmp.py, colrad.py
# 
# The code solves the colisional radiative problem see Stuart loch's 
# spectroscopy class notes or Curt Johnson's upcoming thesis for a derivation
# of the CR matrix.
#
# The read_adf04.py will read adf04 files and produces the dictionary with data
# Takes a density array, temperature array, metastable number, 
#
#
#
################################################################################

import numpy as np
#import fortranformat as ff
import re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
sys.path.append('./')
from r8yip import *
from r8necip import *
from read_adf04_py3 import *


def convert_to_air(lam):
    s = 10**3/lam
    return 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)

def colradpy(fil,metas, temperature_grid, electron_den, use_ionization_in_cr=True, suppliment_with_ecip=True,scale_ecip=-1,use_recombination_three_body= False,use_recombination  = False,td_t=np.array([]),td_n0=np.array([]),err=None,n_trials = 100):
    nsigma=len(metas)
    if(type(fil) == str):
        dict = read_adf04(fil)
    else:
        dict = fil
    if(err):
        err_dict = read_adf04(err)
    
    
    mapp = np.array(['S','P','D','F'])
    nsigmaplus = len(dict['ion_pot'])
    #def colradpy(dict,temperature_grid,electron_den,nsigma=1,
    #             nsigmaplus=1,use_ionization_in_cr=False,suppliment_with_ecip=False)
    #This will calculated the ECIP data for all of the levels and use
    #The ecip where there is no ionization data available in the
    #adf04 file coming from the "S 1   +1" part

    S_ij = np.zeros((len(dict['energy']),len(dict['ion_pot']),len(temperature_grid) ))#changed 4/19/17 to allow mutliple ion metas
    ################################################################################

    if(len(dict['ion_excit']) >0 and use_ionization_in_cr):

        #removed kind='cubic' office computer was throwing error
        #no other computers were throwing errors
        dict['ion_excit'] = np.asarray(dict['ion_excit'])

        ion_excit_interp = interp1d(dict['temp_grid']/11604.5,dict['ion_excit'],axis=1)
        ion_excit_interp_grid = ion_excit_interp(temperature_grid)
        for i in range(0,len(dict['ion_transitions'])):
            for j in range(0,len(temperature_grid)):

                S_ij[dict['ion_transitions'][i,0] -1,dict['ion_transitions'][i,1] -1, j] = ion_excit_interp_grid[i,j]/np.exp(
                    (dict['ion_pot'][dict['ion_transitions'][i,1]-1] - dict['energy'][dict['ion_transitions'][i,0] -1 ])
                    *0.00012398774011749576/temperature_grid[j])


    ################################################################################
    
    if(suppliment_with_ecip):
        dict['ecip'] = np.zeros((len(dict['energy']),len(dict['ion_pot']),len(temperature_grid)))
        for p in range(0,len(dict['ion_pot'])):
            ion_inds = np.where(dict['zpla'][:,p] > -1)[0]#there was a stupid problem beacuse {X} only designates ecip not other ionization
            #ion_inds = np.sort(np.append(ion_inds,dict['ion_transitions'][:,0][np.in1d(np.unique(dict['ion_transitions'][:,0]),ion_inds,invert=True)]))
            ecip_inds = np.where(S_ij[:,p,0] ==0)[0]
            to_use = np.intersect1d(ecip_inds,ion_inds)
            
            dict['ecip'][ion_inds,p,:] =r8ecip(dict['charge_state'], dict['ion_pot'][p],
                             dict['energy'][ion_inds],dict['zpla'][ion_inds,p],temperature_grid*11604.5)

        
            #added this ECIP scaling 12/14/17
            if(scale_ecip>-1):
                # going to do the interpolation between the given data points first

                ion_data_inds = np.sort(dict['ion_transitions'][:,0]) - 1
                                        #-1 for python starting
                                        # at zero
                ion_scale_grid = interp1d(ion_data_inds,S_ij[ion_data_inds]
                                          /dict['ecip'][ion_data_inds],axis=0)

                scaled_ion_interp = ion_scale_grid(np.linspace(ion_data_inds[0],
                           ion_data_inds[-1],ion_data_inds[-1] - ion_data_inds[0] + 1 ))
                dict['ecip_scaled'] = np.zeros_like(dict['ecip'])


                last_ind_scale = ion_data_inds[np.where(scale_ecip > ion_data_inds)[0][-1]]
                dict['ecip_scaled'][0:scale_ecip] = dict['ecip'][0:scale_ecip] *\
                                                    scaled_ion_interp[0:scale_ecip,:]


                dict['ecip_scaled'][scale_ecip:] = dict['ecip'][scale_ecip:]*scaled_ion_interp[scale_ecip,:]

                                                   #+dict['ecip_scaled'][scale_ecip-1] changed 12/18/17
                #now that the interpolation is done go through and do linear extrapolation
                #not sure how valid this is
                ecip_to_use = dict['ecip_scaled'][:,p]
            else:
                ecip_to_use = dict['ecip'][to_use,p,:]
            #for i in range(0,len(to_use)):

            S_ij[to_use,p,:] = ecip_to_use#changed this 12/14/17
                                    #to allow for the scaled ECIP data to be used
                                    #was ecip_to_use[to_use[i],:]


    dict['ionization'] = S_ij
    del S_ij
    
    if(use_recombination):

        recomb_excit_interp = interp1d(dict['temp_grid']/11604.5,
                                       dict['recomb_excit'],
                                       kind='cubic',axis=1)
        
        recomb_excit_interp_grid = recomb_excit_interp(temperature_grid)
        dict['recomb_excit_interp_grid'] = recomb_excit_interp_grid        
        #replace values lower than 1e-30 with a linear interpolation
        #because cubic gives the the wrong answer for some reason

        a,b,c = np.unique(np.where(dict['recomb_excit']<1.e-30)[0],return_inverse=True,return_counts=True)
        #if( temperature_grid < dict['temp_grid'][d.max()]):
        #tmp_interp_grid = []
        if(c.any()):
            if(any(temperature_grid < dict['temp_grid'][c.max()]/11604.5)):
                for v in range(0,len(a)):

                    tmp = np.where(temperature_grid < dict['temp_grid'][c.max()]/11604.5)[0]
                    w = interp1d(dict['temp_grid'][0:c[v]+1]/11604.5,dict['recomb_excit'][a[v],0:c[v]+1],kind='linear')
                    dict['recomb_excit_interp_grid'][a[v],0:c[v]+1] = w(temperature_grid[0:tmp+1])
            
        
            

        dict['recombination'] = np.zeros((len(dict['energy']),len(dict['ion_pot']),len(temperature_grid) ))
        for q in range(0,len(dict['recomb_transitions'])):
            dict['recombination'][dict['recomb_transitions'][q,1]-1,dict['recomb_transitions'][q,0]-1,:] = dict['recomb_excit_interp_grid'][q]

        
    if(use_recombination_three_body):
        #from the los alimos lab manual
        #he['w'] = (he['S']*(he['L']*2+1)-1)/2
        dict['recomb_three_body']=np.zeros((len(dict['energy']),len(dict['ion_pot']),len(temperature_grid)))
        for p in range(0,len(dict['ion_pot'])):
            l_map = np.array(['S','P','D','F','G'])
            if(len(dict['ion_term'][p]) ==2): # was -1)/2.
                w_ion = ( int(dict['ion_term'][p][0]) * (np.where(l_map==dict['ion_term'][p][1])[0][0]*2+1))
            elif(len(dict['ion_term'][p]) ==3):
                w_ion = float(dict['ion_term'][p][2])
            else:
                w_ion=1e30
            dict['w_ion'] = w_ion#1.656742E-22
            dict['recomb_three_body'][:,p,:] = 1.656742E-22* (dict['w'].reshape(len(dict['w']),1)*2+1) /w_ion*np.exp( (dict['ion_pot'][p] - dict['energy']).reshape(len(dict['energy']),1) / (temperature_grid/0.000123985)) *dict['ionization'][:,p,:] / (temperature_grid)**(1.5)
        #rates of this are in cm3
        
    #removed [0:14] 7/7/17 not sure whey that was there
    #removed kind='cubic' 7/17/17, causing singular matrix for adf04-10apr17
    #when run on the office computer no other computer had this error

    
    col_excit_interp = interp1d(dict['temp_grid']/11604.5,dict['col_excit']
                                ,axis=1)






    
    conver = 1.4388
    #indexing of q_ji[j,i] j is the upper level, i is the lower level
    #indexeing of q_ij[i,j], j is the upper level, i is the lower level
    if(err):
        q_ji = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid),n_trials ))
        q_ij = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid),n_trials ))
        A_ji = np.zeros((len(dict['energy']),len(dict['energy']),n_trials))
    else:
        q_ji = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid)))
        q_ij = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid)))
        A_ji = np.zeros((len(dict['energy']),len(dict['energy'])))
        

    if(any(dict['inf_engy'])):
        dict['burg_tully'] = {}
        
        dict['burg_tully']['interp_temp_inds'] = np.where( temperature_grid*11604.5 < dict['temp_grid'][-1])[0]
        dict['burg_tully']['extrap_temp_inds']  = np.where( temperature_grid*11604.5 > dict['temp_grid'][-1])[0]


        FBIG = 0.01
        FZERO = 1E-4
        ELU = np.abs(dict['energy'][dict['col_transitions'][:,0]-1] - dict['energy'][dict['col_transitions'][:,1]-1])/109737.26
        WTU = 2*dict['w'][dict['col_transitions'][:,0]-1]+1
        WTL = 2*dict['w'][dict['col_transitions'][:,1]-1]+1

        S = 3.73491E-10*dict['a_val']*WTU/ELU**3
        FIN = 1/3.*ELU*S/WTL

        
        dict['burg_tully']['type1_ind_arr'] = np.where( (FIN>FBIG) &(dict['S'][dict['col_transitions'][:,0]-1] == dict['S'][dict['col_transitions'][:,1]-1]) & (np.abs(dict['L'][dict['col_transitions'][:,0]-1] - dict['L'][dict['col_transitions'][:,1]-1]) <=1))[0]
        
        dict['burg_tully']['type2_ind_arr'] = np.where( (dict['S'][dict['col_transitions'][:,0]-1] == dict['S'][dict['col_transitions'][:,1]-1]) & ( (np.abs(dict['L'][dict['col_transitions'][:,0]-1] - dict['L'][dict['col_transitions'][:,1]-1]) >1) | (FIN<0.01)))[0]
        

        dict['burg_tully']['type4_ind_arr'] = np.where( (FIN>FZERO) & (dict['S'][dict['col_transitions'][:,0]-1] != dict['S'][dict['col_transitions'][:,1]-1]) & (FIN<FBIG))[0]

        dict['burg_tully']['type3_ind_arr'] = np.where( ((FIN>0.01) | (FIN<FZERO)) & (dict['S'][dict['col_transitions'][:,0]-1] != dict['S'][dict['col_transitions'][:,1]-1]) )[0]

        dict['burg_tully']['c'] = 1.5

        
        def type1_xconv(tconv_grid):
            return 1- np.log(dict['burg_tully']['c']) / np.transpose(np.log( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],1] -1]) + dict['burg_tully']['c']))

        #dict['col_excit'][:,-1] changed 3/11/17   .reshape(len(dict['temp_grid']),1)
        def type1_yconv(tconv_grid):
            return dict['col_excit'][dict['burg_tully']['type1_ind_arr'],:] / np.transpose(np.log( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],1] -1])  + np.exp(1)))

        
        def type2_xconv(tconv_grid):
            return np.transpose(tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type2_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type2_ind_arr'],1] -1])/( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type2_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type2_ind_arr'],1] -1]) +dict['burg_tully']['c']))

        
        def type2_yconv(tconv_grid):
            return dict['col_excit'][dict['burg_tully']['type2_ind_arr'],:]

        
        def type3_xconv(tconv_grid):
            return np.transpose((tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],1] -1]))/( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],1] -1]) +dict['burg_tully']['c']))
        
        def type3_yconv(tconv_grid):
            return np.transpose(tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],1] -1]) +1)*dict['col_excit'][dict['burg_tully']['type3_ind_arr'],:]

        
        def type4_xconv(tconv_grid):
            return 1- np.log(dict['burg_tully']['c']) / np.transpose(np.log( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],1] -1]) + dict['burg_tully']['c']))
        def type4_yconv(tconv_grid):
            return dict['col_excit'][dict['burg_tully']['type4_ind_arr'],:] / np.transpose(np.log( tconv_grid.reshape(len(tconv_grid),1)*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],1] -1])  + dict['burg_tully']['c']))
        


        dict['burg_tully']['type1_xval_arr'] =type1_xconv(dict['temp_grid'])
        dict['burg_tully']['type1_yval_arr'] = type1_yconv(dict['temp_grid'])
        dict['burg_tully']['type1_yval_inf'] = np.abs(dict['inf_engy'][dict['burg_tully']['type1_ind_arr']])

        dict['burg_tully']['type1_coeff_b1'] = (np.log(dict['burg_tully']['type1_yval_inf']) - np.log(dict['burg_tully']['type1_yval_arr'][:,-1])) / (1 - dict['burg_tully']['type1_xval_arr'][:,-1])

        dict['burg_tully']['type1_coeff_a1'] = np.log(dict['burg_tully']['type1_yval_arr'][:,-1]) - dict['burg_tully']['type1_coeff_b1']*dict['burg_tully']['type1_xval_arr'][:,-1]


        
        dict['burg_tully']['type1_coeff_b'] = (np.log(dict['burg_tully']['type1_yval_arr'][:,-1]) - \
                                            np.log(dict['burg_tully']['type1_yval_inf']))/\
                                              (dict['burg_tully']['type1_xval_arr'][:,-1] -1)
        dict['burg_tully']['type1_coeff_a'] = np.exp(np.log(dict['burg_tully']['type1_yval_inf']) -1*dict['burg_tully']['type1_coeff_b'])
        
        type1_zero_inds = np.where(np.isnan(dict['burg_tully']['type1_coeff_a']))[0]
        
        dict['burg_tully']['type1_zero_inds'] = type1_zero_inds

        dict['burg_tully']['type1_coeff_m_lin_fit'] = (np.log(dict['burg_tully']['type1_yval_arr'][dict['burg_tully']['type1_zero_inds'],np.shape(dict['burg_tully']['type1_yval_arr'])[1]-1]) - \
                                      np.log(dict['burg_tully']['type1_yval_arr'][dict['burg_tully']['type1_zero_inds'],-2]))/\
                                       (dict['burg_tully']['type1_xval_arr'][dict['burg_tully']['type1_zero_inds'],np.shape(dict['burg_tully']['type1_yval_arr'])[1]-1] - dict['burg_tully']['type1_xval_arr'][dict['burg_tully']['type1_zero_inds'],-2])
        
        dict['burg_tully']['type1_coeff_b_lin_fit'] = np.log(dict['burg_tully']['type1_yval_arr'][dict['burg_tully']['type1_zero_inds'],-1]) - dict['burg_tully']['type1_coeff_m_lin_fit']*dict['burg_tully']['type1_xval_arr'][dict['burg_tully']['type1_zero_inds'],-1]


        
        if(any(dict['burg_tully']['extrap_temp_inds'])):

            dict['burg_tully']['type1_xval_extrap'] = type1_xconv(temperature_grid[dict['burg_tully']['extrap_temp_inds']]*11604.5)

            #dict['burg_tully']['type1_yval_extrap'] = np.transpose(dict['burg_tully']['type1_coeff_a']*np.exp(np.transpose(dict['burg_tully']['type1_xval_extrap'])*dict['burg_tully']['type1_coeff_b']))
            dict['burg_tully']['type1_yval_extrap'] = np.transpose(np.exp(dict['burg_tully']['type1_coeff_a1'] + np.transpose(dict['burg_tully']['type1_xval_extrap'])*dict['burg_tully']['type1_coeff_b1']))
            
            dict['burg_tully']['type1_excit'] = dict['burg_tully']['type1_yval_extrap']*np.transpose(np.log( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'],1] -1])  + np.exp(1)))



            if(any(dict['burg_tully']['type1_zero_inds'])):

                dict['burg_tully']['type1_yval_extrap_lin'] =  np.transpose(np.exp(dict['burg_tully']['type1_coeff_b_lin_fit'] + dict['burg_tully']['type1_coeff_m_lin_fit']*np.transpose(dict['burg_tully']['type1_xval_extrap'][dict['burg_tully']['type1_zero_inds']])))                

                dict['burg_tully']['type1_excit_lin'] = dict['burg_tully']['type1_yval_extrap_lin']*np.transpose(np.log( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'][dict['burg_tully']['type1_zero_inds']],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type1_ind_arr'][dict['burg_tully']['type1_zero_inds']],1] -1])  + np.exp(1)))

               
        dict['burg_tully']['type2_xval_arr'] =type2_xconv(dict['temp_grid'])
        dict['burg_tully']['type2_yval_arr'] = type2_yconv(dict['temp_grid'])
        dict['burg_tully']['type2_yval_inf'] = np.abs(dict['inf_engy'][dict['burg_tully']['type2_ind_arr']])
        dict['burg_tully']['type2_coeff_b'] = (np.log(dict['burg_tully']['type2_yval_arr'][:,-1]) - \
                                            np.log(dict['burg_tully']['type2_yval_inf']))/\
                                              (dict['burg_tully']['type2_xval_arr'][:,-1] -1)


        dict['burg_tully']['type2_coeff_b1'] = (np.log(dict['burg_tully']['type2_yval_inf']) - np.log(dict['burg_tully']['type2_yval_arr'][:,-1])) / (1 - dict['burg_tully']['type2_xval_arr'][:,-1])

        dict['burg_tully']['type2_coeff_a1'] = np.log(dict['burg_tully']['type2_yval_arr'][:,-1]) - dict['burg_tully']['type2_coeff_b1']*dict['burg_tully']['type2_xval_arr'][:,-1]
        
        dict['burg_tully']['type2_coeff_a'] = np.exp(np.log(dict['burg_tully']['type2_yval_inf']) -1*dict['burg_tully']['type2_coeff_b'])
        dict['burg_tully']['type2_zero_inds'] = np.where(np.isnan(dict['burg_tully']['type2_coeff_a']))[0]
        
        dict['burg_tully']['type2_coeff_m_lin_fit'] = (np.log(dict['burg_tully']['type2_yval_arr'][dict['burg_tully']['type2_zero_inds'],np.shape(dict['burg_tully']['type2_yval_arr'])[1]-1]) - \
                                      np.log(dict['burg_tully']['type2_yval_arr'][dict['burg_tully']['type2_zero_inds'],-2]))/\
                                       (dict['burg_tully']['type2_xval_arr'][dict['burg_tully']['type2_zero_inds'],np.shape(dict['burg_tully']['type2_yval_arr'])[1]-1] - dict['burg_tully']['type2_xval_arr'][dict['burg_tully']['type2_zero_inds'],-2])
        
        dict['burg_tully']['type2_coeff_b_lin_fit'] = np.log(dict['burg_tully']['type2_yval_arr'][dict['burg_tully']['type2_zero_inds'],-1]) - dict['burg_tully']['type2_coeff_m_lin_fit']*dict['burg_tully']['type2_xval_arr'][dict['burg_tully']['type2_zero_inds'],-1]


        if(any(dict['burg_tully']['extrap_temp_inds'])):
            dict['burg_tully']['type2_xval_extrap'] = type2_xconv(temperature_grid[dict['burg_tully']['extrap_temp_inds']]*11604.5)

            if(any(dict['burg_tully']['type2_zero_inds'])):


                dict['burg_tully']['type2_yval_extrap_lin'] =  np.transpose(np.exp(dict['burg_tully']['type2_coeff_b_lin_fit'] + dict['burg_tully']['type2_coeff_m_lin_fit']*np.transpose(dict['burg_tully']['type2_xval_extrap'][dict['burg_tully']['type2_zero_inds']])))                

                dict['burg_tully']['type2_excit_lin'] = dict['burg_tully']['type2_yval_extrap_lin']


            
            #dict['burg_tully']['type2_yval_extrap'] = np.transpose(dict['burg_tully']['type2_coeff_a']*np.exp(np.transpose(dict['burg_tully']['type2_xval_extrap'])*dict['burg_tully']['type2_coeff_b']))
            dict['burg_tully']['type2_yval_extrap'] = np.transpose(np.exp(dict['burg_tully']['type2_coeff_a1'] + np.transpose(dict['burg_tully']['type2_xval_extrap'])*dict['burg_tully']['type2_coeff_b1']))            
            dict['burg_tully']['type2_excit'] = dict['burg_tully']['type2_yval_extrap']

        dict['burg_tully']['type3_xval_arr'] =type3_xconv(dict['temp_grid'])
        dict['burg_tully']['type3_yval_arr'] = type3_yconv(dict['temp_grid'])
        dict['burg_tully']['type3_yval_inf'] = np.abs(dict['inf_engy'][dict['burg_tully']['type3_ind_arr']])




        dict['burg_tully']['type3_coeff_b1'] = (np.log(dict['burg_tully']['type3_yval_inf']) - np.log(dict['burg_tully']['type3_yval_arr'][:,-1])) / (1 - dict['burg_tully']['type3_xval_arr'][:,-1])

        dict['burg_tully']['type3_coeff_a1'] = np.log(dict['burg_tully']['type3_yval_arr'][:,-1]) - dict['burg_tully']['type3_coeff_b1']*dict['burg_tully']['type3_xval_arr'][:,-1]

        
        dict['burg_tully']['type3_coeff_b'] = (np.log(dict['burg_tully']['type3_yval_arr'][:,-1]) - \
                                            np.log(dict['burg_tully']['type3_yval_inf']))/\
                                              (dict['burg_tully']['type3_xval_arr'][:,-1] -1)
        dict['burg_tully']['type3_coeff_a'] = np.exp(np.log(dict['burg_tully']['type3_yval_inf']) -1*dict['burg_tully']['type3_coeff_b'])


        dict['burg_tully']['type3_zero_inds'] = np.where(np.isnan(dict['burg_tully']['type3_coeff_a']))[0]

        dict['burg_tully']['type3_coeff_m_lin_fit'] = (np.log(dict['burg_tully']['type3_yval_arr'][dict['burg_tully']['type3_zero_inds'],np.shape(dict['burg_tully']['type3_yval_arr'])[1]-1]) - \
                                      np.log(dict['burg_tully']['type3_yval_arr'][dict['burg_tully']['type3_zero_inds'],-2]))/\
                                       (dict['burg_tully']['type3_xval_arr'][dict['burg_tully']['type3_zero_inds'],np.shape(dict['burg_tully']['type3_yval_arr'])[1]-1] - dict['burg_tully']['type3_xval_arr'][dict['burg_tully']['type3_zero_inds'],-2])
        
        dict['burg_tully']['type3_coeff_b_lin_fit'] = np.log(dict['burg_tully']['type3_yval_arr'][dict['burg_tully']['type3_zero_inds'],-1]) - dict['burg_tully']['type3_coeff_m_lin_fit']*dict['burg_tully']['type3_xval_arr'][dict['burg_tully']['type3_zero_inds'],-1]


        if(any(dict['burg_tully']['extrap_temp_inds'])):
            dict['burg_tully']['type3_xval_extrap'] = type3_xconv(temperature_grid[dict['burg_tully']['extrap_temp_inds']]*11604.5)
            #dict['burg_tully']['type3_yval_extrap'] = np.transpose(dict['burg_tully']['type3_coeff_a']*np.exp(np.transpose(dict['burg_tully']['type3_xval_extrap'])*dict['burg_tully']['type3_coeff_b']))

            dict['burg_tully']['type3_yval_extrap'] = np.transpose(np.exp(dict['burg_tully']['type3_coeff_a1'] + np.transpose(dict['burg_tully']['type3_xval_extrap'])*dict['burg_tully']['type3_coeff_b1']))
            
            dict['burg_tully']['type3_excit'] = dict['burg_tully']['type3_yval_extrap']/ np.transpose( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'],1] -1])  + 1 )

            if(any(dict['burg_tully']['type3_zero_inds'])):
                
                dict['burg_tully']['type3_yval_extrap_lin'] =  np.transpose(np.exp(dict['burg_tully']['type3_coeff_b_lin_fit'] + dict['burg_tully']['type3_coeff_m_lin_fit']*np.transpose(dict['burg_tully']['type3_xval_extrap'][dict['burg_tully']['type3_zero_inds']])))                

                

                dict['burg_tully']['type3_excit_lin'] = dict['burg_tully']['type3_yval_extrap_lin']/np.transpose( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'][dict['burg_tully']['type3_zero_inds']],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'][dict['burg_tully']['type3_zero_inds']],1] -1])  + 1 )


                #*np.transpose(np.log( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'][dict['burg_tully']['type3_zero_inds']],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type3_ind_arr'][dict['burg_tully']['type3_zero_inds']],1] -1])  + np.exp(1)))


        dict['burg_tully']['type4_xval_arr'] =type4_xconv(dict['temp_grid'])
        dict['burg_tully']['type4_yval_arr'] = type4_yconv(dict['temp_grid'])
        dict['burg_tully']['type4_yval_inf'] = np.abs(dict['inf_engy'][dict['burg_tully']['type4_ind_arr']])
        dict['burg_tully']['type4_coeff_b1'] = (np.log(dict['burg_tully']['type4_yval_inf']) - np.log(dict['burg_tully']['type4_yval_arr'][:,-1])) / (1 - dict['burg_tully']['type4_xval_arr'][:,-1])

        dict['burg_tully']['type4_coeff_a1'] = np.log(dict['burg_tully']['type4_yval_arr'][:,-1]) - dict['burg_tully']['type4_coeff_b1']*dict['burg_tully']['type4_xval_arr'][:,-1]
        
        dict['burg_tully']['type4_coeff_b'] = (np.log(dict['burg_tully']['type4_yval_arr'][:,-1]) - \
                                            np.log(dict['burg_tully']['type4_yval_inf']))/\
                                              (dict['burg_tully']['type4_xval_arr'][:,-1] -1)
        dict['burg_tully']['type4_coeff_a'] = np.exp(np.log(dict['burg_tully']['type4_yval_inf']) -1*dict['burg_tully']['type4_coeff_b'])

        dict['burg_tully']['type4_zero_inds'] = np.where(np.isnan(dict['burg_tully']['type4_coeff_a']))[0]

        dict['burg_tully']['type4_coeff_m_lin_fit'] = (np.log(dict['burg_tully']['type4_yval_arr'][dict['burg_tully']['type4_zero_inds'],np.shape(dict['burg_tully']['type4_yval_arr'])[1]-1]) - \
                                      np.log(dict['burg_tully']['type4_yval_arr'][dict['burg_tully']['type4_zero_inds'],-2]))/\
                                       (dict['burg_tully']['type4_xval_arr'][dict['burg_tully']['type4_zero_inds'],np.shape(dict['burg_tully']['type4_yval_arr'])[1]-1] - dict['burg_tully']['type4_xval_arr'][dict['burg_tully']['type4_zero_inds'],-2])
        
        dict['burg_tully']['type4_coeff_b_lin_fit'] = np.log(dict['burg_tully']['type4_yval_arr'][dict['burg_tully']['type4_zero_inds'],-1]) - dict['burg_tully']['type4_coeff_m_lin_fit']*dict['burg_tully']['type4_xval_arr'][dict['burg_tully']['type4_zero_inds'],-1]

        
        if(any(dict['burg_tully']['extrap_temp_inds'])):
            dict['burg_tully']['type4_xval_extrap'] = type4_xconv(temperature_grid[dict['burg_tully']['extrap_temp_inds']]*11604.5)
            #dict['burg_tully']['type4_yval_extrap'] = np.transpose(dict['burg_tully']['type4_coeff_a']*np.exp(np.transpose(dict['burg_tully']['type4_xval_extrap'])*dict['burg_tully']['type4_coeff_b']))
            dict['burg_tully']['type4_yval_extrap'] = np.transpose(np.exp(dict['burg_tully']['type4_coeff_a1'] + np.transpose(dict['burg_tully']['type4_xval_extrap'])*dict['burg_tully']['type4_coeff_b1']))            
            dict['burg_tully']['type4_excit'] = dict['burg_tully']['type4_yval_extrap']*np.transpose(np.log( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'],1] -1])  + dict['burg_tully']['c']))
            if(any(dict['burg_tully']['type4_zero_inds'])):
                
                dict['burg_tully']['type4_yval_extrap_lin'] =  np.transpose(np.exp(dict['burg_tully']['type4_coeff_b_lin_fit'] + dict['burg_tully']['type4_coeff_m_lin_fit']*np.transpose(dict['burg_tully']['type4_xval_extrap'][dict['burg_tully']['type4_zero_inds']])))

                dict['burg_tully']['type4_excit_lin'] = dict['burg_tully']['type4_yval_extrap_lin']*np.transpose(np.log( temperature_grid[dict['burg_tully']['extrap_temp_inds']].reshape(len(dict['burg_tully']['extrap_temp_inds']),1)*11604.5*0.69488623/ (dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'][dict['burg_tully']['type4_zero_inds']],0]-1] - dict['energy'][dict['col_transitions'][dict['burg_tully']['type4_ind_arr'][dict['burg_tully']['type4_zero_inds']],1] -1])  + dict['burg_tully']['c']))


        #col_excit_interp_grid_tmp = 
        col_excit_interp_grid = np.zeros((len(dict['col_excit']),len(temperature_grid)))

        for i in range(0,len(dict['burg_tully']['extrap_temp_inds'])):
            col_excit_interp_grid[dict['burg_tully']['type1_ind_arr'],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type1_excit'][:,i]
            col_excit_interp_grid[dict['burg_tully']['type2_ind_arr'],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type2_excit'][:,i]

            col_excit_interp_grid[dict['burg_tully']['type3_ind_arr'],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type3_excit'][:,i]

                
            col_excit_interp_grid[dict['burg_tully']['type4_ind_arr'],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type4_excit'][:,i]

            if(any(dict['burg_tully']['type1_zero_inds'])):
                col_excit_interp_grid[dict['burg_tully']['type1_ind_arr'][dict['burg_tully']['type1_zero_inds']],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type1_excit_lin'][:,i]
            if(any(dict['burg_tully']['type2_zero_inds'])):
                col_excit_interp_grid[dict['burg_tully']['type2_ind_arr'][dict['burg_tully']['type2_zero_inds']],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type2_excit_lin'][:,i]
            if(any(dict['burg_tully']['type3_zero_inds'])):
                col_excit_interp_grid[dict['burg_tully']['type3_ind_arr'][dict['burg_tully']['type3_zero_inds']],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type3_excit_lin'][:,i]
            if(any(dict['burg_tully']['type4_zero_inds'])):
                col_excit_interp_grid[dict['burg_tully']['type4_ind_arr'][dict['burg_tully']['type4_zero_inds']],dict['burg_tully']['extrap_temp_inds'][i] ] = dict['burg_tully']['type4_excit_lin'][:,i]
                
        col_excit_interp_grid[:,dict['burg_tully']['interp_temp_inds']] = col_excit_interp(temperature_grid[dict['burg_tully']['interp_temp_inds']])

        
    else:
        col_excit_interp_grid = col_excit_interp(temperature_grid)



    dict['col_excit_interp'] = col_excit_interp_grid

    #adding in error prop stuff
    if(err):
        err_excit_interp = interp1d(dict['temp_grid']/11604.5,err_dict['col_excit']
                                            ,axis=1)
        dict['err_excit_interp'] = err_excit_interp(temperature_grid)
        dict['err_excit_trials'] = np.zeros((np.shape(dict['err_excit_interp'])[0],np.shape(dict['err_excit_interp'])[1],n_trials))
        for i in range(0,len(temperature_grid)):
                                              
            dict['err_excit_trials'][:,i,:] = np.transpose(np.random.multivariate_normal(dict['col_excit_interp'][:,0],np.diag(dict['err_excit_interp'][:,0]),size=n_trials))
    
    rates_ji = np.zeros_like(col_excit_interp_grid)
    rates_ij = np.zeros_like(col_excit_interp_grid)
    
    for i in range(0,len(dict['col_transitions'])):
        # the if statement here is to account for files that have had energies
        # deleted from the begining of the adf04. It is assumed if those
        # entries are not at the top of the adf04 file then the user does not
        # want to included then this allows the matrix to be made ignoring
        # levels above the levels taken out in the file.
        if( dict['col_transitions'][i,0] < len(dict['energy'])+1):
            
            for j in range(0,len(temperature_grid)):
                if(err):
                    q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j,:] = 1/(2*dict['w'][dict['col_transitions'][i,0]-1] +1)* np.sqrt(13.6058/temperature_grid[j]) * dict['err_excit_trials'][i,j]*2.1716E-8
                    
                    q_ij[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j,:] = (2*dict['w'][dict['col_transitions'][i,0]-1]+1)/(2*dict['w'][dict['col_transitions'][i,1]-1]+1) * np.exp(-np.abs(dict['energy'][dict['col_transitions'][i,0]-1] -dict['energy'][dict['col_transitions'][i,1]-1])/8.065E3/temperature_grid[j]) *q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,j,:]
                    A_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,:] = dict['a_val'][i]
                    
                else:
                    q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j] = 1/(2*dict['w'][dict['col_transitions'][i,0]-1] +1)* np.sqrt(13.6058/temperature_grid[j]) * col_excit_interp_grid[i,j]*2.1716E-8
                    rates_ji[i,j] = 1/(2*dict['w'][dict['col_transitions'][i,0]-1] +1)* np.sqrt(13.6058/temperature_grid[j]) * col_excit_interp_grid[i,j]*2.1716E-8

                    q_ij[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j] = (2*dict['w'][dict['col_transitions'][i,0]-1]+1)/(2*dict['w'][dict['col_transitions'][i,1]-1]+1) * np.exp(-np.abs(dict['energy'][dict['col_transitions'][i,0]-1] -dict['energy'][dict['col_transitions'][i,1]-1])/8.065E3/temperature_grid[j]) *q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,j]
                
                    rates_ij[i,j] = (2*dict['w'][dict['col_transitions'][i,0]-1]+1)/(2*dict['w'][dict['col_transitions'][i,1]-1]+1) * np.exp(-np.abs(dict['energy'][dict['col_transitions'][i,0]-1] -dict['energy'][dict['col_transitions'][i,1]-1])/8.065E3/temperature_grid[j]) *q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,j]

                    A_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1] = dict['a_val'][i]

    #transpose the q_ij matrix so the indexes make sense
    if(err):
        q_ij = q_ij.transpose(1,0,2,3)        
    else:
        q_ij = q_ij.transpose(1,0,2)

    #construct the collisional radiative matrix
    if(use_recombination):
        cr = np.zeros( (len(dict['energy']),len(dict['energy'])+nsigmaplus,len(temperature_grid),len(electron_den)  )  )

    else:
        if(err):
            cr = np.zeros( (len(dict['energy']),len(dict['energy']),len(temperature_grid),len(electron_den) ,n_trials )  )
        else:
            cr = np.zeros( (len(dict['energy']),len(dict['energy']),len(temperature_grid),len(electron_den)  )  )

    cr_loss = np.zeros_like(cr)
    
    for t in range(0,len(temperature_grid)):
        for e in range(0,len(electron_den)):
            for i in range(0,len(dict['energy'])):
                for mm in range(0,n_trials):
                    #level i depopulating mechanisms
                    #these are all the transitions from the level
                    cr[i,i,t,e,mm]  =  -1*np.sum(A_ji[i,:])
                    cr_loss[i,i,t,e,mm]  =  -1*np.sum(A_ji[i,:])
                    #these are all the excitation and dexcitation bringing you out of the level
                    cr[i,i,t,e,mm] = cr[i,i,t,e,mm] - np.sum(q_ji[i,:,t,mm])*electron_den[e] - np.sum(q_ij[i,:,t,mm])*electron_den[e]
                    cr_loss[i,i,t,e,mm] = cr_loss[i,i,t,e,mm] - np.sum(q_ji[i,:,t,mm])*electron_den[e] - np.sum(q_ij[i,:,t,mm])*electron_den[e]                
                #these are the ways to ionize out of ion
                    if(use_ionization_in_cr):
                        cr[i,i,t,e,mm] = cr[i,i,t,e,mm] - np.sum(dict['ionization'][i,:,t])*electron_den[e]
                        cr_loss[i,i,t,e,mm] = cr_loss[i,i,t,e,mm] - np.sum(dict['ionization'][i,:,t])*electron_den[e]                    

                    #
                    #level i populating mechanisms
                    #
                    #these are the transition rates from higher levels into the level i
                    cr[i,0:len(dict['energy']),t,e,mm] = cr[i,0:len(dict['energy']),t,e,mm] + A_ji[:,i]
                    #these are excitation and dexciation into the level i
                    cr[i,0:len(dict['energy']),t,e,mm] = cr[i,:len(dict['energy']),t,e,mm] + q_ij[:,i,t,mm]*electron_den[e] + q_ji[:,i,t,mm]*electron_den[e]

                if(use_recombination):
                    qq=0
                    for p in range(0,len(dict['ion_pot'])):

                        cr[0:len(dict['energy']),len(dict['energy'])+p,t,e] = cr[0:len(dict['energy']),len(dict['energy'])+p,t,e] + dict['recombination'][0:len(dict['energy']),p,t]*electron_den[e]#recomb_excit_interp_grid[0:len(dict['energy']),t][qq:(p+1)*len(dict['energy'])]*electron_den[e]
                    #qq = (p+1)*len(dict['energy'])-1
                if(use_recombination_three_body):
                    for p in range(0,len(dict['ion_pot'])):
                        cr[0:len(dict['energy']),len(dict['energy'])+p,t,e] = cr[0:len(dict['energy']),len(dict['energy'])+p,t,e] + dict['recomb_three_body'][0:len(dict['energy']),p,t]*(electron_den[e])**2

            #these are the ionizations into the
            #cr[-1,0:len(dict['energy']),t,e] = cr[-1,0:len(dict['energy']),t,e] + dict['ionization'][:,t]*electron_den[e]

    #now its time for the reduced matrix
    #recombination and would go here
    #also octopole transitions
    #del q_ij
    #del q_ji
    dict['col_rates_ij'] = rates_ij
    dict['col_rates_ji'] = rates_ji
    dict['cr'] = cr
    levels_to_keep = np.setdiff1d( np.linspace(0,len(dict['energy'])-1,
                                               len(dict['energy']) ,dtype='int64'),metas)
    if(use_recombination):
        beta_tmp = -cr[nsigma:len(dict['energy']),metas]
        
        #beta_tmp = np.append(-cr[nsigma:len(dict['energy']),metas],#was 0:sigma
        #                     -cr[nsigma:len(dict['energy']),len(dict['energy']),:,:][:,None,:,:],axis=1)
        #                    #was -cr[nsigma:len(dict['energy']),0,:,:][:,None,:,:],axis=1)

        for p in range(0,nsigmaplus):
            beta_tmp = np.append(beta_tmp,-cr[nsigma:len(dict['energy']),len(dict['energy'])+p,:,:][:,None,:,:],axis=1)
            
        aa_tmp = cr[nsigma:len(dict['energy']),nsigma:len(dict['energy'])]
        aa_tmp_inv = np.zeros((len(aa_tmp),len(aa_tmp), len(temperature_grid),len(electron_den)))
    else:
        beta_tmp = -cr[levels_to_keep][:,metas]
        aa_tmp = cr[levels_to_keep][:,levels_to_keep]
        aa_tmp_inv = np.zeros((len(aa_tmp),len(aa_tmp), len(temperature_grid),len(electron_den)))
        if(td_t.any() and td_n0.any()):
            dict['td_n0'] = td_n0
            dict['td_t'] = td_t
            if(use_ionization_in_cr):
                dict['td_cr'] = np.zeros( (len(dict['energy'])+1,len(dict['energy'])+1,len(temperature_grid),len(electron_den)  )  )

                dict['td_cr'][-1,0:len(dict['energy']),:,:] = np.einsum('ij,k->ijk',dict['ionization'][:,0,:],electron_den)
            else:
                dict['td_cr'] = np.zeros( (len(dict['energy']),len(dict['energy']),len(temperature_grid),len(electron_den)  )  )
                
            dict['td_cr'][0:len(dict['energy']),0:len(dict['energy']),:,:] = dict['cr']

            dict['td_pop'] = np.zeros( (len(dict['td_n0']),len(dict['td_t']),len(temperature_grid),len(electron_den)))
            for e in range(0,len(electron_den)):
                for t in range(0,len(temperature_grid)):
                    eigenval,eigenvectors = np.linalg.eig(dict['td_cr'][:,:,t,e])
                    v0 = np.dot(np.linalg.inv(eigenvectors),dict['td_n0'])

                    vt = v0[:,None]*np.exp(eigenval[:,None]*dict['td_t'])
                    dict['td_pop'][:,:,t,e] = np.dot(eigenvectors,vt)
            dict['eigenval'] = eigenval
            dict['eigenvectors'] = eigenvectors





    if(err):
        aa_tmp_inv = np.linalg.inv(aa_tmp.transpose(2,3,4,0,1)).transpose(3,4,0,1,2)
        aa_inv = aa_tmp_inv

    else:
        for i in range(0,len(electron_den)):
            for j in range(0,len(temperature_grid)):

                aa_tmp_inv[:,:,j,i] = np.linalg.inv(aa_tmp[:,:,j,i])
                aa_inv = aa_tmp_inv


    if(use_recombination):
        populations = np.zeros((len(aa_tmp),nsigma+nsigmaplus,len(temperature_grid),len(electron_den) ))
    else:
        if(err):
            populations = np.zeros((len(aa_tmp),nsigma,len(temperature_grid),len(electron_den) ,n_trials))            
        else:
            populations = np.zeros((len(aa_tmp),nsigma,len(temperature_grid),len(electron_den) ))
    for t in range(0,len(temperature_grid)):
        for e in range(0,len(electron_den)):
            for n in range(0,n_trials):
                populations[:,:,t,e,n] = np.dot(aa_inv[:,:,t,e,n],beta_tmp[:,:,t,e,n])
    dict['populations'] = populations

    pecs = []
    sbxs = []
    pecs_levels = []
    wavelengths =[]
    specific_line_pwr = []
    if(nsigma>1):
        driving_population_norm=False
    else:
        driving_population_norm=True
    for i in range(0,len(A_ji)):
        for j in range(nsigma,len(A_ji)):

            if( levels_to_keep[j-nsigma]>i and A_ji[j,i].any() >1E-31):

                if(driving_population_norm):

                    pecs.append(A_ji[j,i]*populations[j-nsigma]/electron_den/
                                (1+np.sum(populations,axis=0)))
                    specific_line_pwr.append(A_ji[j,i]*populations[j-nsigma]/electron_den/
                                            (1+np.sum(populations,axis=0))*((dict['energy'][j] - dict['energy'][i])/5.03e15))                    
                else:
                    pecs.append(A_ji[j,i]*populations[j-nsigma]/electron_den)
                    specific_line_pwr.append(A_ji[j,i]*populations[j-nsigma]/electron_den\
                                            *((dict['energy'][j] - dict['energy'][i])/5.03e15))                    
                    
                pecs_levels.append(np.array([levels_to_keep[j-nsigma],i]))

                wavelengths.append(  (1./abs(dict['energy'][levels_to_keep[j-nsigma]] - dict['energy'][i])*1e7))#/
                                     #convert_to_air(1./abs(dict['energy'][j] - dict['energy'][i])*1e7))
                #rad_pwr = rad_pwr + (he['energy'][j] - he['energy'][i])*(A_ji[j,i]*populations[j-nsigma])
                #print(i,j)

    wavelengths = np.asarray(wavelengths)
    pecs = np.asarray(pecs)
    dict['pecs'] = pecs
    pecs_levels = np.asarray(pecs_levels)



    ################################################################################
    #
    #Branching ratios
    #meta stable cross coupling 
    #
    ################################################################################

    qcd = np.zeros((nsigma,nsigma-1,len(temperature_grid),len(electron_den),n_trials))

    poptmp = np.zeros((len(aa_tmp),len(aa_tmp),nsigma,len(temperature_grid),
                       len(electron_den),n_trials ))
    poptmpr = np.zeros((len(aa_tmp),len(aa_tmp),nsigmaplus,len(temperature_grid),
                        len(electron_den),n_trials ))

    scd = np.zeros((nsigma,nsigmaplus,len(temperature_grid),len(electron_den),n_trials))
    acd = np.zeros((nsigma,nsigmaplus,len(temperature_grid),len(electron_den),n_trials))
    xcd = np.zeros((nsigmaplus,nsigmaplus-1,len(temperature_grid),len(electron_den),n_trials))

    for j in range(0,nsigma):
        for k in range(0,len(temperature_grid)):
            for l in range(0,len(electron_den)):        
                for i in range(0,len(beta_tmp)):
                    for mm in range(0,n_trials):
                        
                        poptmp[:,i,j,k,l,mm] = aa_inv[:,i,k,l,mm]*beta_tmp[i,j,k,l,mm]

                        #if there is only one metastable we have to normalize to the entire population
                        #this is because there can be more populations in the excited states than in the
                        #ground state. This was only implimented after 8/29/17 as the problem was noticed
                        #looking at wlike_mons11#w0.dat with only one metstable.
                        if(driving_population_norm):
                            poptmp[:,:,j,k,l,mm] = poptmp[:,:,j,k,l,mm]/(1+ np.sum(np.sum(poptmp[:,:,j,k,l,mm],axis=1),axis=0)) 

                        F = np.sum(poptmp[:,:,j,k,l,mm],axis=1)
                        #if(use_recombination):
                            #R = np.sum(aa_inv[:,i,k,l]*dict['recombination'][i,k])#recomb_excit_interp_grid[i,k])

                        for m in range(0,nsigmaplus):
                            if(use_ionization_in_cr):

                                #The commented section is wrong it is now fixed and checked for he is 1 and 2 metastables
                                #not sure why I was doing /(1 + np.sum(np.sum(    this doesn't make sense
                                #for some reason it was impacting the metastable case but not the ground only... weird
                                #scd[j,m,k,l] = dict['ionization'][j,k]/(1+ np.sum(np.sum(poptmp[:,:,j,k,l]))) +np.sum(
                                #    dict['ionization'][nsigma:,k]*F)
                                if( driving_population_norm):

                                    scd[j,m,k,l,mm] = dict['ionization'][j,m,k]/(1+np.sum(populations[:,j,k,l,mm]))+np.sum(dict['ionization'][nsigma:,m,k]*F)
                                else:
                                    scd[j,m,k,l],mm = dict['ionization'][metas[j],m,k]+np.sum(dict['ionization'][levels_to_keep,m,k]*F)
                            if(use_recombination):

                                poptmpr[:,i,m,k,l,mm] =aa_inv[i,:,k,l,mm]*beta_tmp[i,nsigma+m,k,l,mm]/electron_den[l]
                                #acd[j,m,k,l] = (dict['recomb_three_body'][j,m,k]*electron_den[l] + dict['recombination'][j,m,k]) /(1+np.sum(populations[:,j,k,l])) + np.sum(cr[j,nsigma:len(dict['energy']),k,l]*np.sum(aa_inv[:,:,k,l]*dict['recombination'][nsigma:,m,k]*dict['recomb_three_body'][j,m,k]*electron_den[l],axis=1)) #(dict['recomb_three_body'][j,k]*electron_den[l] + dict['recomb_excit_interp_grid'][k,j]) /(1+np.sum(populations[:,j,k,l])) - np.sum(F * (dict['recomb_three_body'][nsigma:,k]*electron_den[l] + dict['recomb_excit_interp_grid'][nsigma:,k
                                acd[j,m,k,l,mm] = (dict['recomb_three_body'][j,m,k,mm]*electron_den[l] + dict['recombination'][j,m,k,mm])  - np.sum(
                                    cr[j,nsigma:len(dict['energy']),k,l,mm]*np.sum(aa_inv[:,:,k,l,mm]*(dict['recomb_three_body'][nsigma:,m,k,mm]*electron_den[l] + dict['recombination'][nsigma:,m,k,mm]),axis=1))

                                mindd=0
                                for n in range(0,nsigmaplus-1):
                                    if(m !=n ):
                                        xcd[m,mindd,k,l,mm] = -np.sum(dict['ionization'][nsigma:,n,k,mm]*np.sum(aa_inv[:,:,k,l,mm]*(dict['recomb_three_body'][nsigma:,m,k,mm]*electron_den[l] + dict['recombination'][nsigma:,m,k,mm]),axis=1))
                                        mindd = mindd+1

                        mind=0
                        #for n in range(0,nsigma):
                        for n in metas:

                            if(n !=metas[j]):
                                qcd[j,mind,k,l,mm] = (cr[n,metas[j],k,l,mm] + np.sum(cr[n,levels_to_keep,k,l,mm]*F))/electron_den[l]
                                #cr_del = np.delete(cr,metas,axis=1)
                                #qcd[j,mind,k,l] = (cr[n,j,k,l] + np.sum(cr_del[n,:,k,l]*F))/electron_den[l]
                                mind = mind+1


            '''
            poptmprec = np.zeros((len(aa_tmp),len(aa_tmp),nsigma,nsigmaplus,len(temperature_grid),
                               len(electron_den) ))
            '''

        '''
        if( use_recombination):
            for j in range(0,nsigma):
                for m in range(0,nsigmaplus):
                    for k in range(0,len(temperature_grid)):
                        for l in range(0,len(electron_den)):        
                            for i in range(0,len(beta_tmp)):
                                poptmprec[:,i,j,m,k,l] = aa_inv[:,i,k,l]*beta_tmp[i,nsigma+m,k,l]/electron_den[l]
                            R = np.sum(poptmprec[:,:,j,m,k,l],axis=1)

                            acd[j,m,k,l] = (dict['recombination'][j,m,k]+dict['recomb_three_body'][j,m,k]*electron_den[l]) - np.sum(cr[n,nsigma:len(dict['energy']),k,l]*R)



            dict['poptmprec'] = poptmprec

        '''




        #scdd = 0
        #for ii in range(0,len(aa_inv)):
        #        for jj in range(0,len(aa_inv)):
        #                scdd = scdd + dict['ionization'][i]*np.sum(aa_inv[j,i,0,0]*beta_tmp[i,0,0,0])
        #scdd = scdd + dict['ionization'][0,0]
        #print(scdd)
    
    
    sxbs = np.zeros_like(pecs)
    for m in range(0,nsigma):
        for t in range(0,len(temperature_grid)):
            for n in range(0,len(electron_den)):
                for mm in range(0,n_trials):
                    sxbs[:,m,t,n,mm] = scd[m,0,t,n,mm]/pecs[:,m,t,n,mm]

    dict['aa_inv'] = aa_inv
    dict['aa_tmp'] = aa_tmp
    dict['beta_tmp'] = beta_tmp
    dict['a_ji'] = A_ji
    dict['populations'] = populations
    dict['sxbs'] = sxbs
    dict['acd'] = acd
    dict['scd'] = scd
    dict['qcd'] = qcd
    dict['xcd'] = xcd
    dict['poptmp'] = poptmp
    dict['poptmpr'] = poptmpr
    dict['F'] = F
    dict['nsigma'] = nsigma
    dict['user_temp_grid'] = temperature_grid
    dict['user_dens_grid'] = electron_den
    dict['pecs'] = pecs
    dict['specific_line_pwr'] = np.asarray(specific_line_pwr)
    dict['wavelengths'] = wavelengths
    dict['loss'] = cr_loss
    dict['pecs_levels'] = pecs_levels
    dict['plt_pwr'] = np.sum(dict['specific_line_pwr'],axis=0)
    dict['metas'] = np.asarray(metas)
    return dict


#def normalize_ecip(wdx):
#    ion_data_inds = np.sort(wdx['ion_transitions'][:,0] - 1#-1 for python starting
                                                           # at zero




'''

    def convert_adas_config_to_nist_config(dict):
        #first part remove all of the trailing 1's NIST does not use them
        connection_levels = pymysql.connect(host='localhost',
                                            user='root',
                                            password='spectra',
                                            db='levels',
                                            charset='utf8mb4',
                                            cursorclass=pymysql.cursors.DictCursor)
        cur_levels = connection_levels.cursor()
        cur_levels.execute(
           'SELECT * FROM ASD_Levels WHERE element= %s AND spectr_charge=%s',
            (dict['element'], dict['roman_numeral']))# get all of the energy levels
                                               #for a given element and charge state
        b = cur_levels.fetchone()# b is now the dictionary with all the data from NIST
        nist_shells = re.split('[.]',str(b['conf']))

        

        adas_shells =  re.split('[.]',wdx['config'][0])
        index_adas_start = 0
        for i in range(0,len(adas_shells)):
            if(adas_shells[i] == nist_shells[0]):
                index_adas_start = i
        trim = index_adas_start
        self.nist_config = self.config
        print(adas_shells,nist_shells)
        for i in range(0,len(self.config)):
            split_arr = re.split('[.]',self.config[i])
            print(split_arr)
            for j in range(0,len(split_arr)):
                if(split_arr[j][-1] == '1'):
                    split_arr[j] = split_arr[j][0:len(split_arr[j]) -1]

            tmp_config_str = split_arr[trim]
            for k in range(trim+1,len(split_arr)):
                tmp_config_str = tmp_config_str + '.' + split_arr[k]
            self.nist_config[i] = tmp_config_str





connection_levels = pymysql.connect(host='localhost',
                                     user='root',
                                     password='spectra',
                                     db='levels',
                                     charset='utf8mb4',
                                     cursorclass=pymysql.cursors.DictCursor)
cur_levels = connection_levels.cursor()
cur_levels.execute(
   'SELECT * FROM ASD_Levels WHERE element= %s AND spectr_charge=%s',
    (self.element, self.roman_numeral))# get all of the energy levels
                                               #for a given element and charge state
b = cur_levels.fetchall()

for i in range(0,len(wdx1['config'])):
    r = len(wdx1['config'][i][4:])/3
    tmp=''
    for j in range(0,r):


trim=1
confs = []
for i in range(0,len(wdx1['config'])):
    split_arr = re.split('[ ]',wdx1['config'][i])
    print(split_arr)
    for j in range(0,len(split_arr)):
        if(split_arr[j][-1] == '1'):
            split_arr[j] = split_arr[j][0:len(split_arr[j]) -1]

        tmp_config_str = split_arr[trim]
        for k in range(trim+1,len(split_arr)):
            tmp_config_str = tmp_config_str + '.' + split_arr[k]
        confs.append(tmp_config_str)

'''

'''

for i in range(0,4):
        for j in range(0,4):
                plt.figure()
                plt.xlabel('level')
                plt.ylabel('ionization percent from level')
                plt.scatter(np.linspace(1,250,249),c['ionization'][1:,i]*np.sum(c['poptmp'][:,:,0,i,j],axis=1)/c['scd'][0,0,i,j])
                plt.title('Temperature: ' + str(c['user_temp_grid'][i]) + '    Density: ' + str(c['user_dens_grid'][j]) )


'''








'''
if(dict['charge_state'] == 0):
ion_bal = np.zeros((nsigma+1,nsigma+1))
ion_bal[0,0] = -1*(qcd[0,0,0,0]+scd[0,0,0,0]) 
ion_bal[1,1] = -1*( qcd[1,0,0,0] +scd[1,0,0,0])
ion_bal[1,0] = qcd[0,0,0,0]
ion_bal[0,1] = qcd[1,0,0,0]

ion_bal[0,2] = acd[0,0,0,0]
ion_bal[1,2] = acd[1,0,0,0]
ion_bal[2,2] = -1*np.sum(acd[:,0,0,0])
ion_bal[2,0] = scd[0,0,0,0]
ion_bal[2,1] = scd[1,0,0,0]
'''



def plot_lev(poptmp,leve,meta,temp,dens):
    plt.scatter(np.linspace(0,len(poptmp)-1,len(poptmp)),poptmp[leve,:,meta,temp,dens]/np.sum(poptmp[leve,:,meta,temp,dens]),label=str(i),color=np.random.rand(3,1))
    plt.legend()
    plt.title(str(leve))
    plt.xlabel('level')
    plt.ylabel('Populating percent')
'''
def stuart_plot(lev,norm):
q = []
co = ['b','g','r','k','c','y']
for i in range(0,6):
    q.append(np.sum(poptmp[lev,:,i,:,0],axis=0))
fig = plt.figure()
for i in range(0,6):
    plt.plot(temperature_grid,q[i],color=co[i],label=i)
plt.figure()
plt.plot(np.linspace(0,len(poptmp),len(poptmp)-1), poptmp[lev,:,norm,0,0])


pec1 = 397
pec2 = 419
temp = 7
den = 3

def pec_dec(pec1,pec2,temp,den):
fig = plt.figure()
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)

ax1.plot(electron_den,np.sum(pecs[pec1,:,temp,:],axis=0)/np.sum(pecs[pec2,:,temp,:],axis=0))
ax1.set_xscale('log')
ax2.plot(temperature_grid,np.sum(pecs[pec1,:,:,den],axis=0)/np.sum(pecs[pec2,:,:,den],axis=0))
metas_cont = np.sum(poptmp[pecs_levels[pec1,0]-nsigma,:,:,temp,den],axis=0)
pop_ways = []
for i in range(0,nsigma):
    pop_ways.append(metas_cont[i]*poptmp[pecs_levels[pec1,0]-nsigma,:,i,temp,den])
    ax3.scatter(np.linspace(0,len(dict['energy']) -1-nsigma,len(dict['energy'])-nsigma),pop_ways[i])
ax3.set_ylim(0,np.max(pop_ways)*1.1)
metas_cont2 = np.sum(poptmp[pecs_levels[pec2,0]-nsigma,:,:,temp,den],axis=0)    
pop_ways2 = []
for i in range(0,nsigma):
    pop_ways2.append(metas_cont2[i]*poptmp[pecs_levels[pec2,0]-nsigma,:,i,temp,den])
    ax4.scatter(np.linspace(0,len(dict['energy']) -1-nsigma,len(dict['energy'])-nsigma),pop_ways2[i])
ax4.set_ylim(0,np.max(pop_ways2)*1.2)

for i in range(0,nsigma):
    ax5.plot(temperature_grid,np.sum(poptmp[pecs_levels[pec1,0]-nsigma,:,i,:,den],axis=0))
    ax6.plot(temperature_grid,np.sum(poptmp[pecs_levels[pec2,0]-nsigma,:,i,:,den],axis=0))
ax5.set_yscale('log')
ax6.set_yscale('log')





In [90]: wdxm['pecs'][541,0,0,0]/(np.sum(wdx['populations'])+1) + np.sum(wdxm['pecs'][541,1:6,0,0]/(np.sum(wdx['populations'])+1)*wdx['populations'][0:5,0,0,0])


In [95]: wdxm['scd'][0,0,0,0]/(np.sum(wdx['populations'])+1) + np.sum(wdxm['scd'][1:6,0,0,0]/(np.sum(wdx['populations'])+1)*wdx['populations'][0:5,0,0,0])

'''





