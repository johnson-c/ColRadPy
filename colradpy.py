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
# todo     allow non consequative metastable levels, remove read adf04 part
#          at the begining and replace with read_adf04.py, put into a class
#
#
################################################################################

import numpy as np
import fortranformat as ff
import re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
sys.path.append('./')
from r8yip import *
from r8necip import *
from read_adf04 import *
nsigmaplus = 1

def convert_to_air(lam):
    s = 10**3/lam
    return 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)

def colradpy(fil,metas, temperature_grid, electron_den, use_ionization_in_cr=True, suppliment_with_ecip=True,scale_ecip=-1,use_recombination_three_body= True,use_recombination  = True):
    nsigma=len(metas)
    if(type(fil) == str):
        dict = read_adf04(fil)
    else:
        dict = fil
    mapp = np.array(['S','P','D','F'])
    #def colradpy(dict,temperature_grid,electron_den,nsigma=1,
    #             nsigmaplus=1,use_ionization_in_cr=False,suppliment_with_ecip=False)
    #This will calculated the ECIP data for all of the levels and use
    #The ecip where there is no ionization data available in the
    #adf04 file coming from the "S 1   +1" part

    S_ij = np.zeros((len(dict['energy']),len(temperature_grid) ))
    ################################################################################

    if(len(dict['ion_excit']) >0 and use_ionization_in_cr):

        #removed kind='cubic' office computer was throwing error
        #no other computers were throwing errors
        dict['ion_excit'] = np.asarray(dict['ion_excit'])
        ion_excit_interp = interp1d(dict['temp_grid']/11604.5,dict['ion_excit'],axis=1)
        ion_excit_interp_grid = ion_excit_interp(temperature_grid)
        for i in range(0,len(dict['ion_transitions'])):
            for j in range(0,len(temperature_grid)):

                S_ij[dict['ion_transitions'][i,0] -1,j] = ion_excit_interp_grid[i,j]/np.exp(
                    (dict['ion_pot'] - dict['energy'][dict['ion_transitions'][i,0] -1 ])
                    *0.00012398774011749576/temperature_grid[j])

    ################################################################################

    if(suppliment_with_ecip):

        ion_inds = np.where(dict['zpla'] > -1)[0]
        dict['ecip'] =r8ecip(dict['charge_state'], dict['ion_pot'],
                             dict['energy'][ion_inds],dict['zpla'][ion_inds],temperature_grid*11604.5)
        
        ecip_inds = np.where(S_ij[:,0] ==0)[0]
        to_use = np.intersect1d(ecip_inds,ion_inds)
        
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
            ecip_to_use = dict['ecip_scaled']
        else:
            ecip_to_use = dict['ecip']
        
        for i in range(0,len(to_use)):
            S_ij[to_use[i],:] = ecip_to_use[to_use[i],:]#changed this 12/14/17
                                    #to allow for the scaled ECIP data to be used

    dict['ionization'] = S_ij
    del S_ij
    
    if(use_recombination):

        recomb_excit_interp = interp1d(dict['temp_grid']/11604.5,
                                       dict['recomb_excit'],
                                       kind='cubic',axis=1)
        recomb_excit_interp_grid = recomb_excit_interp(temperature_grid)
        dict['recomb_excit_interp_grid'] = recomb_excit_interp_grid
    if(use_recombination_three_body):
        #from the los alimos lab manual
        #he['w'] = (he['S']*(he['L']*2+1)-1)/2
        l_map = np.array(['S','P','D','F','G'])
        if(len(dict['ion_term']) ==2): # was -1)/2.
            w_ion = ( int(dict['ion_term'][0]) * (np.where(l_map==dict['ion_term'][1])[0][0]*2+1))
        elif(len(dict['ion_term']) ==3):
            w_ion = float(dict['ion_term'][2])
        else:
            w_ion=1e30
        dict['w_ion'] = w_ion#1.656742E-22
        dict['recomb_three_body'] = 1.656742E-22* (dict['w'].reshape(len(dict['w']),1)*2+1) /w_ion*np.exp( (dict['ion_pot'] - dict['energy']).reshape(len(dict['energy']),1) / (temperature_grid/0.000123985)) *dict['ionization'] / (temperature_grid)**(1.5)
        #rates of this are in cm3
        
    #removed [0:14] 7/7/17 not sure whey that was there
    #removed kind='cubic' 7/17/17, causing singular matrix for adf04-10apr17
    #when run on the office computer no other computer had this error
    col_excit_interp = interp1d(dict['temp_grid']/11604.5,dict['col_excit']
                                ,axis=1)

    conver = 1.4388
    #indexing of q_ji[j,i] j is the upper level, i is the lower level
    #indexeing of q_ij[i,j], j is the upper level, i is the lower level
    q_ji = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid) ))
    q_ij = np.zeros((len(dict['energy']),len(dict['energy']),len(temperature_grid) ))
    A_ji = np.zeros((len(dict['energy']),len(dict['energy'])))
    col_excit_interp_grid = col_excit_interp(temperature_grid)
    dict['col_excit_interp'] = col_excit_interp_grid
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

                q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j] = 1/(2*dict['w'][dict['col_transitions'][i,0]-1] +1)* np.sqrt(13.6058/temperature_grid[j]) * col_excit_interp_grid[i,j]*2.1716E-8
                rates_ji[i,j] = 1/(2*dict['w'][dict['col_transitions'][i,0]-1] +1)* np.sqrt(13.6058/temperature_grid[j]) * col_excit_interp_grid[i,j]*2.1716E-8

                q_ij[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1, j] = (2*dict['w'][dict['col_transitions'][i,0]-1]+1)/(2*dict['w'][dict['col_transitions'][i,1]-1]+1) * np.exp(-np.abs(dict['energy'][dict['col_transitions'][i,0]-1] -dict['energy'][dict['col_transitions'][i,1]-1])/8.065E3/temperature_grid[j]) *q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,j]
                
                rates_ij[i,j] = (2*dict['w'][dict['col_transitions'][i,0]-1]+1)/(2*dict['w'][dict['col_transitions'][i,1]-1]+1) * np.exp(-np.abs(dict['energy'][dict['col_transitions'][i,0]-1] -dict['energy'][dict['col_transitions'][i,1]-1])/8.065E3/temperature_grid[j]) *q_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1,j]

            A_ji[dict['col_transitions'][i,0]-1,dict['col_transitions'][i,1]-1] = dict['a_val'][i]

    #transpose the q_ij matrix so the indexes make sense
    q_ij = q_ij.transpose(1,0,2)

    #construct the collisional radiative matrix
    if(use_recombination):
        cr = np.zeros( (len(dict['energy']),len(dict['energy'])+1,len(temperature_grid),len(electron_den)  )  )

    else:
        cr = np.zeros( (len(dict['energy']),len(dict['energy']),len(temperature_grid),len(electron_den)  )  )

    cr_loss = np.zeros_like(cr)
    
    for t in range(0,len(temperature_grid)):
        for e in range(0,len(electron_den)):
            for i in range(0,len(dict['energy'])):
                #level i depopulating mechanisms
                #these are all the transitions from the level
                cr[i,i,t,e]  =  -1*np.sum(A_ji[i,:])
                cr_loss[i,i,t,e]  =  -1*np.sum(A_ji[i,:])
                #these are all the excitation and dexcitation bringing you out of the level
                cr[i,i,t,e] = cr[i,i,t,e] - np.sum(q_ji[i,:,t])*electron_den[e] - np.sum(q_ij[i,:,t])*electron_den[e]
                cr_loss[i,i,t,e] = cr_loss[i,i,t,e] - np.sum(q_ji[i,:,t])*electron_den[e] - np.sum(q_ij[i,:,t])*electron_den[e]                
                #these are the ways to ionize out of ion
                if(use_ionization_in_cr):
                    cr[i,i,t,e] = cr[i,i,t,e] - dict['ionization'][i,t]*electron_den[e]
                    cr_loss[i,i,t,e] = cr_loss[i,i,t,e] - dict['ionization'][i,t]*electron_den[e]                    

                #
                #level i populating mechanisms
                #
                #these are the transition rates from higher levels into the level i
                cr[i,0:len(dict['energy']),t,e] = cr[i,0:len(dict['energy']),t,e] + A_ji[:,i]
                #these are excitation and dexciation into the level i
                cr[i,0:len(dict['energy']),t,e] = cr[i,:len(dict['energy']),t,e] + q_ij[:,i,t]*electron_den[e] + q_ji[:,i,t]*electron_den[e]

            if(use_recombination):
                cr[0:len(dict['energy']),-1,t,e] = cr[0:len(dict['energy']),-1,t,e] + recomb_excit_interp_grid[0:len(dict['energy']),t]*electron_den[e]
            if(use_recombination_three_body):

                cr[0:len(dict['energy']),-1,t,e] = cr[0:len(dict['energy']),-1,t,e] + dict['recomb_three_body'][0:len(dict['energy']),t]*(electron_den[e])**2

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
        beta_tmp = np.append(-cr[nsigma:len(dict['energy']),0:nsigma],
                             -cr[nsigma:len(dict['energy']),-1,:,:][:,None,:,:],axis=1)
        aa_tmp = cr[nsigma:len(dict['energy']),nsigma:len(dict['energy'])]
        aa_tmp_inv = np.zeros((len(aa_tmp),len(aa_tmp), len(temperature_grid),len(electron_den)))

    else:

        beta_tmp = -cr[levels_to_keep][:,metas]
        aa_tmp = cr[levels_to_keep][:,levels_to_keep]
        aa_tmp_inv = np.zeros((len(aa_tmp),len(aa_tmp), len(temperature_grid),len(electron_den)))

    for i in range(0,len(electron_den)):
        for j in range(0,len(temperature_grid)):

            aa_tmp_inv[:,:,j,i] = np.linalg.inv(aa_tmp[:,:,j,i])
            aa_inv = aa_tmp_inv


    if(use_recombination):
        populations = np.zeros((len(aa_tmp),nsigma+1,len(temperature_grid),len(electron_den) ))
    else:
        populations = np.zeros((len(aa_tmp),nsigma,len(temperature_grid),len(electron_den) ))
    for t in range(0,len(temperature_grid)):
        for e in range(0,len(electron_den)):
            populations[:,:,t,e] = np.dot(aa_inv[:,:,t,e],beta_tmp[:,:,t,e])


    pecs = []
    sbxs = []
    pecs_levels = []
    wavelengths =[]
    rad_pwr = np.zeros( (nsigma,len(temperature_grid),len(electron_den) ))
    if(nsigma>1):
        driving_population_norm=False
    else:
        driving_population_norm=True


    for i in range(0,len(A_ji)):
        for j in range(nsigma,len(A_ji)):
            if( j>i and A_ji[j,i] >1E-31):
                if(driving_population_norm):

                    pecs.append(A_ji[j,i]*populations[j-nsigma]/electron_den/
                                (1+np.sum(populations,axis=0)))
                    
                else:
                    pecs.append(A_ji[j,i]*populations[j-nsigma]/electron_den)
                    #print('a')
                #rad_pwr = rad_pwr + (energy[j] - energy[i])*A[j,i]*populations[j-nsigma]
                pecs_levels.append(np.array([j,i]))

                wavelengths.append(  (1./abs(dict['energy'][j] - dict['energy'][i])*1e7))#/
                                     #convert_to_air(1./abs(dict['energy'][j] - dict['energy'][i])*1e7))
                #rad_pwr = rad_pwr + (he['energy'][j] - he['energy'][i])*(A_ji[j,i]*populations[j-nsigma])
                #print(i,j)
                
    wavelengths = np.asarray(wavelengths)
    pecs = np.asarray(pecs)
    pecs_levels = np.asarray(pecs_levels)



    ################################################################################
    #
    #Branching ratios
    #meta stable cross coupling 
    #
    ################################################################################

    qcd = np.zeros((nsigma,nsigma-1,len(temperature_grid),len(electron_den)))

    poptmp = np.zeros((len(aa_tmp),len(aa_tmp),nsigma,len(temperature_grid),
                       len(electron_den) ))

    scd = np.zeros((nsigma,nsigmaplus,len(temperature_grid),len(electron_den)))
    acd = np.zeros((nsigma,nsigmaplus,len(temperature_grid),len(electron_den)))
    for j in range(0,nsigma):
        for k in range(0,len(temperature_grid)):
            for l in range(0,len(electron_den)):        
                for i in range(0,len(beta_tmp)):
                    poptmp[:,i,j,k,l] = aa_inv[:,i,k,l]*beta_tmp[i,j,k,l]

                #if there is only one metastable we have to normalize to the entire population
                #this is because there can be more populations in the excited states than in the
                #ground state. This was only implimented after 8/29/17 as the problem was noticed
                #looking at wlike_mons11#w0.dat with only one metstable.
                if(driving_population_norm):
                    poptmp[:,:,j,k,l] = poptmp[:,:,j,k,l]/(1+ np.sum(np.sum(poptmp[:,:,j,k,l],axis=1),axis=0)) 

                F = np.sum(poptmp[:,:,j,k,l],axis=1)

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
                            scd[j,m,k,l] = dict['ionization'][j,k]/(1+np.sum(populations[:,j,k,l]))+np.sum(dict['ionization'][nsigma:,k]*F)
                        else:
                            scd[j,m,k,l] = dict['ionization'][j,k]+np.sum(dict['ionization'][nsigma:,k]*F)
                    if(use_recombination):

                        acd[j,m,k,l] = (dict['recomb_three_body'][j,k]*electron_den[l] + dict['recomb_excit_interp_grid'][j,k]) /(1+np.sum(populations[:,j,k,l])) + np.sum(cr[j,nsigma:len(dict['energy']),k,l]*np.sum(aa_inv[:,:,k,l]*dict['recomb_excit_interp_grid'][nsigma:,k]*dict['recomb_three_body'][j,k]*electron_den[l],axis=1)) #(dict['recomb_three_body'][j,k]*electron_den[l] + dict['recomb_excit_interp_grid'][k,j]) /(1+np.sum(populations[:,j,k,l])) - np.sum(F * (dict['recomb_three_body'][nsigma:,k]*electron_den[l] + dict['recomb_excit_interp_grid'][nsigma:,k]))
                    '''
                        if(use_recombination_three_body):
                            acd[j,m,k,l] = (recomb_excit_interp_grid[j,k]+dict['recomb_three_body'][j,k]*electron_den[l]) - np.sum(
                                cr[j,nsigma:len(dict['energy']),k,l]*np.sum(
                                    aa_inv[:,:,k,l]*recomb_excit_interp_grid[nsigma:,k]*dict['recomb_three_body'][j,k]*electron_den[l],
                                                     axis=1))
                        else:    
                            acd[j,m,k,l] = recomb_excit_interp_grid[j,k] - np.sum(
                                cr[j,nsigma:len(dict['energy']),k,l]*np.sum(
                                    aa_inv[:,:,k,l]*recomb_excit_interp_grid[nsigma:,k],
                                    axis=1))

                    else:
                        if(use_recombination_three_body):

                            acd[j,m,k,l] = (dict['recomb_three_body'][j,k]*electron_den[l]) - np.sum(
                                cr[j,nsigma:len(dict['energy']),k,l]*np.sum(
                                    aa_inv[:,:,k,l]*dict['recomb_three_body'][j,k]*electron_den[l],
                                                     axis=1))

                '''
                mind=0
                for n in range(0,nsigma):
                    if(n !=j):
                        qcd[j,mind,k,l] = (cr[n,j,k,l] + np.sum(cr[n,nsigma:len(dict['energy']),k,l]*F))/electron_den[l]
                        mind = mind+1


    poptmprec = np.zeros((len(aa_tmp),len(aa_tmp),nsigma,nsigmaplus,len(temperature_grid),
                       len(electron_den) ))

    for j in range(0,nsigma):
        for m in range(0,nsigmaplus):
            for k in range(0,len(temperature_grid)):
                for l in range(0,len(electron_den)):        
                    for i in range(0,len(beta_tmp)):

                        poptmprec[:,i,j,m,k,l] = aa_inv[:,i,k,l]*beta_tmp[i,-1,k,l]
                    R = np.sum(poptmprec[:,:,j,m,k,l],axis=1)
                    acd[j,m,k,l] = ((dict['recomb_three_body'][j,k]*electron_den[l]**2 + dict['recomb_excit_interp_grid'][j,k]*electron_den[l]) + np.sum( cr[j,nsigma:len(dict['energy']),k,l]*R) )/electron_den[l]


    dict['poptmprec'] = poptmprec






                        
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
                sxbs[:,m,t,n] = scd[m,0,t,n]/pecs[:,m,t,n]

    dict['aa_inv'] = aa_inv
    dict['aa_tmp'] = aa_tmp
    dict['beta_tmp'] = beta_tmp
    dict['a_ji'] = A_ji
    dict['populations'] = populations
    dict['sxbs'] = sxbs
    dict['acd'] = acd
    dict['scd'] = scd
    dict['qcd'] = qcd
    dict['poptmp'] = poptmp
    dict['F'] = F
    dict['nsigma'] = nsigma
    dict['user_temp_grid'] = temperature_grid
    dict['user_dens_grid'] = electron_den
    dict['pecs'] = pecs
    dict['wavelengths'] = wavelengths
    dict['loss'] = cr_loss
    dict['pecs_levels'] = pecs_levels
    dict['rad_pwr'] = rad_pwr
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
