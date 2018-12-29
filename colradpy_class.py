################################################################################
# file name         : colrad.py
# author            : Curt Johnson
# description       : This code solves the collsion radiative matrix
# version           : 0.3
# python version    : Python 3.6.5 ipython 6.4.0 anaconda
# dependencies      : read_adf04_py3_class.py, numpy, scipy interpolate,
#                   : matplotlib, re, sys, r8yip, r8necip, ecip_rates
#                   : burgess_tully_rates, collections
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
import re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
sys.path.append('./')
from r8yip import *
from r8necip import *
from read_adf04_py3_class import *
from ecip_rates import *
from burgess_tully_rates import *
from split_multiplet import *
import collections

def convert_to_air(lam):
    """This function converts the vacuum wavelength of spectral lines to 
       the air wavelength. The wavelength is the difference in energy
       levels of the upper and lower transition. By defualt in the adf04
       file these energy correspond to the vacuum wavelengths to this
       must be converted at some point for UV and visible light.

    :param lam: The wavelengths in nm in vacuum
    :type lam: array

    :returns: array --wavelengths in nm in air

    """
    
    s = 10**3/lam
    return 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)

class colradpy():

    """The ColRadPy class, this class will store data and carry out calculation to solve
       the collisional radiative set of equations. A full tutorial is provided in the 
       documentation.
   

    :param fil: the file path to the input file
    :type fil: string

    :param metas: array of the levels that metastable
    :type metas: integer array

    :param temp_grid: array of temperature for calculation (eV)
    :type metas: float array

    :param dens_grid: array of temperature for calculation (cm-3)
    :type metas: float array

    :param use_ionization: Flag to turn on ionization in calculation, default = True
    :type metas: bool

    :param suppliment_with_ecip: Flag to turn on ECIP supplimentation to ionization rates, default = Tue
    :type metas: bool

    :param use_recombination_three_body: Flag to turn 3 body recombination on,  default = True
    :type metas: bool

    :param use_recombination: Flag to turn recombination on, default = True
    :type metas: bool

    :param td_t: Time array for time dependent solution to CR equations
    :type metas: float array

    :param td_n0: Initial populations of levels at time t=0 for TD CR equations
    :type metas: float array

    :param td_source: source term for populations in the TD CR equations
    :type metas: float array

    """
    def __init__(self,fil,metas=np.array([]),temp_grid=np.array([]),electron_den=np.array([]),
                 use_ionization=True,suppliment_with_ecip=True,use_recombination_three_body=True,
                 use_recombination = True, td_t = np.array([]), td_n0=np.array([]),td_source=np.array([]),
                  default_pop_norm=True):

        self.data = {}
        self.processed = {} 
        self.data['user'] = {}
        self.data['user']['temp_grid'] = temp_grid
        self.data['user']['dens_grid'] = electron_den
        self.data['user']['use_ionization'] = use_ionization
        self.data['user']['suppliment_with_ecip'] = suppliment_with_ecip
        self.data['user']['use_recombination_three_body'] = use_recombination_three_body
        self.data['user']['use_recombination'] = use_recombination
        self.data['user']['metas'] = 'look at [atomic][metas] for values'
        self.data['user']['td_t'] = td_t
        self.data['user']['td_n0'] = td_n0
        self.data['user']['td_source'] = td_source
        self.data['user']['default_pop_norm'] = default_pop_norm
        self.populate_data(fil)
        self.data['atomic']['metas'] = np.asarray(metas)

    def update_dict(self,d, u):
        """This function will update dictionary that is stores all parameter for the class

        :param d: the dictionary that is to be updated
        :type d: dictionary

        :param w: the dictionary the use wants to update to
        :type w: dictionary
        """

        for k, v in u.items():
            if isinstance(v, collections.Mapping):
                d[k] = self.update_dict(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    def populate_data(self,fil):
        """This function will populate atomic data. Currently this only uses the
           adf04 file but there is nothing special about this format and should be
           put out to pasture at some point

        :param fil: path to the file that contains atomic data
        :type fil: file path

        :returns: array --wavelengths in nm in air

        """
        if(type(fil) == str):

            self.data = self.update_dict(self.data,read_adf04(fil))
            self.data['user']['file_loc'] = fil
        else:
            if(self.data):
                self.data = self.update_dict(self.data,fil)
            else:
                self.data = fil

    def make_ecip(self):
        """This function calls the ecip_rates function and populates the ecip
           values in the dictionary. See documentation ecip_rates.py for a better desciption
           
        """

        self.data['rates']['ioniz']['ecip'] = np.zeros((len(self.data['atomic']['energy']),
                                                        len(self.data['atomic']['ion_pot']),
                                                        len(self.data['user']['temp_grid'])))
        #there was a stupid problem beacuse {X} only designates ecip not other ionization

        self.data['rates']['ioniz']['ecip'] = ecip_rates(self.data['atomic']['energy'],
                                                  self.data['atomic']['ion_pot'],self.data['atomic']['zpla'],
                                                  self.data['atomic']['zpla1'],self.data['atomic']['charge_state'],
                                                  self.data['user']['temp_grid'])
        
    def make_burgess_tully(self):
        """This function calls the burgess_tully_rates function and updates the 'rates' dictionary

           values in the dictionary. See documentation ecip_rates.py for a better desciption
           
        """
        self.data['rates'].update( burgess_tully_rates(self.data['user']['temp_grid'],self.data['input_file']['temp_grid'],
                                                       self.data['rates']['excit']['col_transitions'],self.data['rates']\
                                                       ['excit']['col_excit'],
                                                       self.data['atomic']['energy'],self.data['atomic']['w'],
                                                       self.data['rates']['a_val'],self.data['atomic']['S'],
                                                       self.data['atomic']['L'],self.data['rates']['inf_engy']))

    def make_ioniz_from_reduced_ionizrates(self):

        """This function calls will make ionization rates to be used in the CR matrix from the
           reduced ionization rates that are provided in the adf04 file. This function
           must be called even if ionization is not provided so that ECIP rates can be 
           supplimented into the matrix
        """

        
        #if there is no ionization in the file create a zero matrix so that ECIP can be subbed in
        if(np.size(self.data['rates']['ioniz']['ion_excit']) < 1):
            print('No ionization in the input file ECIP can be made')
            self.data['rates']['ioniz']['ionization'] = np.zeros((len(self.data['atomic']['energy']),
                                                                  len(self.data['atomic']['ion_pot']),
                                                                  len(self.data['user']['temp_grid'])))
            return

        #changed 4/19/17 to allow mutliple ion metas        
        self.data['rates']['ioniz']['ionization'] = np.zeros((len(self.data['atomic']['energy']),
                                                                       len(self.data['atomic']['ion_pot']),
                                                                       len(self.data['user']['temp_grid'])))

        ion_excit_interp = interp1d(self.data['input_file']['temp_grid']/11604.5,
                                    self.data['rates']['ioniz']['ion_excit'],axis=1,kind='slinear')
        
        ion_excit_interp_grid = ion_excit_interp(self.data['user']['temp_grid'])

        for i in range(0,len(self.data['rates']['ioniz']['ion_transitions'])):
            for j in range(0,len(self.data['user']['temp_grid'])):

                self.data['rates']['ioniz']['ionization'][self.data['rates']['ioniz']\
                    ['ion_transitions'][i,0] -1,self.data['rates']['ioniz']\
                    ['ion_transitions'][i,1] -1, j] = ion_excit_interp_grid[i,j]/\
                    np.exp((self.data['atomic']['ion_pot'][self.data['rates']['ioniz']\
                    ['ion_transitions'][i,1]-1] - self.data['atomic']['energy']\
                    [self.data['rates']['ioniz']['ion_transitions'][i,0] -1 ])\
                    *0.00012398774011749576/self.data['user']['temp_grid'][j])

    def suppliment_with_ecip(self):
        """This function will suppliment the ionization that is to be used in the CR matrix
           with ECIP ionization if there is no ionization that is provided. While ECIP ionization
           is better than not including ionization using other ionization rates is generall better
           Ionization must be must be included to accurately model plasma even if it is approximate ECIP
        """

        if('ecip' not in self.data['rates']['ioniz']):
            print("ECIP was not previously calculated, calculating now"+\
                   " from inside 'suppliment_with_ecip'")
            self.make_ecip()

        if(('ionization' not in self.data['rates']['ioniz']) and (self.data['user']['use_ionization'])):
            print('No ionization rates from files were calculated,'+\
                  ' try to calculate them')
            self.make_ioniz_from_reduced_ionizrates()

        #added for the case that ionization is turned off and 
        if('ionization' not in self.data['rates']['ioniz']):
            self.data['rates']['ioniz']['ionization'] = np.zeros((len(self.data['atomic']['energy']),
                                                                  len(self.data['atomic']['ion_pot']),
                                                                  len(self.data['user']['temp_grid'])))
        
        for p in range(0,len(self.data['atomic']['ion_pot'])):
            #there was a stupid problem beacuse {X} only designates ecip not other ionization   
            ion_inds = np.where( self.data['atomic']['zpla'][:,p] > -1)[0]
            ion_inds2 = np.where( self.data['atomic']['zpla1'] == p +1)[0]
            ion_inds = np.intersect1d(ion_inds,ion_inds2)
            ecip_inds = np.where(self.data['rates']['ioniz']['ionization'][:,p,0] ==0)[0]
            to_use = np.intersect1d(ecip_inds,ion_inds)
            self.data['rates']['ioniz']['ionization'][to_use,p,:] =\
                                        self.data['rates']['ioniz']['ecip'][to_use,p,:]

    def make_recombination_rates_from_file(self):
        """This function will make recombination rates from the rates that are provided in the
           adf04 file.
        """

        recomb_excit_interp = interp1d(self.data['input_file']['temp_grid']/11604.5,
                                       self.data['rates']['recomb']['recomb_excit'],
                                       axis=1,kind='slinear')
        self.data['rates']['recomb']['recomb_excit_interp_grid'] =\
                                    recomb_excit_interp(self.data['user']['temp_grid'])
        #replace values lower than 1e-30 with a linear interpolation
        #because slinear gives the the wrong answer for some reason
        #maybe dont need this part now, 1dinterp was using kind='slinear'
        #changed now and hopefully wont get the errors
        a,b,c = np.unique(np.where(self.data['rates']['recomb']['recomb_excit']\
                                   <1.e-30)[0],return_inverse=True,return_counts=True)
        if(c.any()):
            if(any(self.data['user']['temp_grid'] < \
                   self.data['input_file']['temp_grid'][c.max()]/11604.5)):
                for v in range(0,len(a)):
                    tmp = np.where(self.data['user']['temp_grid'] <\
                                   self.data['input_file']['temp_grid'][c.max()]/11604.5)[0]
                    w = interp1d(self.data['input_file']['temp_grid'][0:c[v]+1]/11604.5,\
                                 self.data['rates']['recomb']['recomb_excit'][a[v],0:c[v]+1],\
                                 kind='slinear')

                    self.data['rates']['recomb']['recomb_excit_interp_grid'][a[v],0:c[v]+1] =\
                                                    w(self.data['user']['temp_grid'][0:tmp+1])

        self.data['rates']['recomb']['recombination'] = \
            np.zeros((len(self.data['atomic']['energy']),
                      len(self.data['atomic']['ion_pot']),
                      len(self.data['user']['temp_grid'])))
        
        for q in range(0,len(self.data['rates']['recomb']['recomb_transitions'])):
            self.data['rates']['recomb']['recombination']\
                [self.data['rates']['recomb']['recomb_transitions'][q,1]-1,
                 self.data['rates']['recomb']['recomb_transitions'][q,0]-1,:] = \
                                    self.data['rates']['recomb']['recomb_excit_interp_grid'][q]

    def make_three_body_recombination(self):
        """This function will make three body recombination rates by using a detailed balance
           relation with ionization that used. Three body recombination becomes important at high densities
           and low temperatures
        """
        
        #from the los alimos lab manual
        #he['w'] = (he['S']*(he['L']*2+1)-1)/2
        self.data['rates']['recomb']['recomb_three_body'] = \
                        np.zeros((len(self.data['atomic']['energy']),
                                  len(self.data['atomic']['ion_pot']),
                                  len(self.data['user']['temp_grid'])))
        for p in range(0,len(self.data['atomic']['ion_pot'])):
            l_map = np.array(['S','P','D','F','G'])
            if(len(self.data['atomic']['ion_term'][p]) ==2): # was -1)/2.
                w_ion = ( int(self.data['atomic']['ion_term'][p][0]) * \
                          (np.where(l_map==self.data['atomic']['ion_term'][p][1])[0][0]*2+1))
            elif(len(self.data['atomic']['ion_term'][p]) ==3):
                w_ion = float(self.data['atomic']['ion_term'][p][2])
            else:
                w_ion=1e30
            self.data['rates']['recomb']['w_ion'] = w_ion#1.656742E-22
            #this is from detailed balance
            self.data['rates']['recomb']['recomb_three_body'][:,p,:] = 1.656742E-22* \
                (self.data['atomic']['w'].reshape(len(self.data['atomic']['w']),1)*2+1)/ \
                w_ion*np.exp( (self.data['atomic']['ion_pot'][p] - \
                self.data['atomic']['energy']).reshape(len(self.data['atomic']['energy']),1) / \
                (self.data['user']['temp_grid']/0.000123985)) *\
                self.data['rates']['ioniz']['ionization'][:,p,:] / \
                (self.data['user']['temp_grid'])**(1.5)
        #rates of this are in cm3

    def make_electron_excitation_rates(self):

        """This function will make both electron impact excitation and deexcitation rates
           from the values that stored in adf04 file on a user defined temperature grid
           If values are above the last calculate point in the adf04 then a burgess tully
           extrapolation will be used. There is currently no extrapolation below the first
           calculated temperature point. THis is something to add in the future.
        """
        tmp = interp1d(self.data['input_file']['temp_grid']/11604.5,self.data['rates']['excit']['col_excit'],axis=1,kind='slinear')

        self.data['rates']['excit']['col_excit_interp'] = np.zeros((len(self.data['rates']['excit']['col_excit']),
                                                                         len(self.data['user']['temp_grid'])))


        if( ('burg_tully' not in self.data['rates'].keys()) and (np.max(self.data['input_file']\
                                ['temp_grid'])/11604.5 < np.max(self.data['user']['temp_grid']) )):
            print('Atleast one user temp pont above last calculated temperature using extrapolation be carefull')
            self.make_burgess_tully()

        if(np.max(self.data['input_file']['temp_grid'])/11604.5 < \
           np.max(self.data['user']['temp_grid'])):

            for i in range(0,len(self.data['rates']['burg_tully']['extrap_temp_inds'])):
                for j in range(0,3):
                    self.data['rates']['excit']['col_excit_interp']\
                        [self.data['rates']['burg_tully']['ind_arrs'][j],
                                      self.data['rates']['burg_tully']['extrap_temp_inds'][i] ] =\
                                      self.data['rates']['burg_tully']['excit_extrap'][j][:,i]

                for j in range(0,3):
                    if(self.data['rates']['burg_tully']['zero_inds'][j].size >0):
                        self.data['rates']['excit']['col_excit_interp']\
                            [self.data['rates']['burg_tully']['ind_arrs'][j]\
                                           [self.data['rates']['burg_tully']['zero_inds'][j]],
                                        self.data['rates']['burg_tully']['extrap_temp_inds'][i] ] =\
                                        self.data['rates']['burg_tully']['excit_extrap_lin'][j][:,i]


            self.data['rates']['excit']['col_excit_interp'][:,self.data['rates']['burg_tully']\
                                   ['interp_temp_inds']] =\
            tmp(self.data['user']['temp_grid'][self.data['rates']['burg_tully']['interp_temp_inds']])

        else:

            self.data['rates']['excit']['col_excit_interp'] = tmp(self.data['user']['temp_grid'])


    def populate_cr_matrix(self):
        """This function will populate the collision radiative matrix with all the rates that
           user asks to be included into the calculation.
        """

        if('col_excit_interp' not in self.data['rates']['excit'].keys()):
            print('Electron collisional rates have not yet been made on the user defined temperature grid. Doing that now')
            self.make_electron_excitation_rates()
        self.data['cr_matrix'] = {}
        self.data['cr_matrix']['q_ji'] = np.zeros((len(self.data['atomic']['energy']),
                                                   len(self.data['atomic']['energy']),
                                                   len(self.data['user']['temp_grid'])))
        self.data['cr_matrix']['q_ij'] = np.zeros((len(self.data['atomic']['energy']),
                                                   len(self.data['atomic']['energy']),
                                                   len(self.data['user']['temp_grid'])))
        self.data['cr_matrix']['A_ji'] = np.zeros((len(self.data['atomic']['energy']),
                                                   len(self.data['atomic']['energy'])))
        for i in range(0,len(self.data['rates']['excit']['col_transitions'])):
            # the if statement here is to account for files that have had energies
            # deleted from the begining of the adf04. It is assumed if those
            # entries are not at the top of the adf04 file then the user does not
            # want to included then this allows the matrix to be made ignoring
            # levels above the levels taken out in the file.
            if( self.data['rates']['excit']['col_transitions'][i,0] < len(self.data['atomic']['energy'])+1):
                for j in range(0,len(self.data['user']['temp_grid'])):

                    self.data['cr_matrix']['q_ji'][self.data['rates']['excit']['col_transitions'][i,0]-1,
                                                   self.data['rates']['excit']['col_transitions'][i,1]-1, j] =  \
                    1/(2*self.data['atomic']['w'][self.data['rates']['excit']['col_transitions'][i,0]-1] +1)* \
                    np.sqrt(13.6058/self.data['user']['temp_grid'][j]) * \
                    self.data['rates']['excit']['col_excit_interp'][i,j]*2.1716E-8

                    self.data['cr_matrix']['q_ij'][self.data['rates']['excit']['col_transitions'][i,0]-1,
                                                   self.data['rates']['excit']['col_transitions'][i,1]-1, j] =  \
                    (2*self.data['atomic']['w'][self.data['rates']['excit']['col_transitions'][i,0]-1]+1)/  \
                    (2*self.data['atomic']['w'][self.data['rates']['excit']['col_transitions'][i,1]-1]+1) * \
                    np.exp(-np.abs(self.data['atomic']['energy'][self.data['rates']['excit']['col_transitions'][i,0]-1] - \
                    self.data['atomic']['energy'][self.data['rates']['excit']['col_transitions'][i,1]-1])/  \
                    8.065E3/self.data['user']['temp_grid'][j]) * \
                    self.data['cr_matrix']['q_ji'][self.data['rates']['excit']['col_transitions'][i,0]-1,
                                                   self.data['rates']['excit']['col_transitions'][i,1]-1,j]

                self.data['cr_matrix']['A_ji'][self.data['rates']['excit']['col_transitions'][i,0]-1,
                                               self.data['rates']['excit']['col_transitions'][i,1]-1] = \
                                                                                            self.data['rates']['a_val'][i]
        #transpose the self.data['cr_matrix']['q_ij'] matrix so the indexes make sense
        self.data['cr_matrix']['q_ij'] = self.data['cr_matrix']['q_ij'].transpose(1,0,2)

        
        nsigmaplus = len(self.data['atomic']['ion_pot'])
        if(self.data['user']['use_recombination']):
            if(len(np.unique(self.data['rates']['recomb']['recomb_transitions'][:,0]))>nsigmaplus):
               nsigmaplus = len(np.unique(self.data['rates']['recomb']['recomb_transitions'][:,0]))
        elif(self.data['user']['use_recombination_three_body']):
            nsigmaplus = np.shape(self.data['rates']['recomb']['recomb_three_body'])[1]
        
        self.data['cr_matrix']['cr'] = np.zeros((len(self.data['atomic']['energy'])+nsigmaplus,
                                                 len(self.data['atomic']['energy'])+nsigmaplus,
                                                 len(self.data['user']['temp_grid']),
                                                 len(self.data['user']['dens_grid'])))

        self.data['cr_matrix']['cr_loss'] = np.zeros((len(self.data['atomic']['energy'])+nsigmaplus,
                                                 len(self.data['atomic']['energy'])+nsigmaplus,
                                                 len(self.data['user']['temp_grid']),
                                                 len(self.data['user']['dens_grid'])))


        
        for i in range(0,len(self.data['atomic']['energy'])):
            #level i depopulating mechanisms
            #these are all the transitions from the level

            self.data['cr_matrix']['cr'][i,i,:,:]  =  -1*np.sum(self.data['cr_matrix']['A_ji'][i,:])
            self.data['cr_matrix']['cr_loss'][i,i,:,:]  =  -1*np.sum(self.data['cr_matrix']['A_ji'][i,:])
            #these are all the excitation and dexcitation bringing you out of the level
            self.data['cr_matrix']['cr'][i,i,:,:] = self.data['cr_matrix']['cr'][i,i,:,:] - \
            np.sum(np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ji'][i,:,:],
                             self.data['user']['dens_grid']),axis=0)-\
            np.sum(np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ij'][i,:,:],
                             self.data['user']['dens_grid']),axis=0)

            self.data['cr_matrix']['cr_loss'][i,i,:,:] = self.data['cr_matrix']['cr'][i,i,:,:] - \
            np.sum(np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ji'][i,:,:],
                             self.data['user']['dens_grid']),axis=0)-\
            np.sum(np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ij'][i,:,:],
                             self.data['user']['dens_grid']),axis=0)
            
            #these are the ways to ionize out of ion
            if(self.data['rates']['ioniz']['ionization'].size > 0):
                self.data['cr_matrix']['cr'][i,i,:,:] = self.data['cr_matrix']['cr'][i,i,:,:] - \
                                      np.sum(np.einsum('ij,k->ijk',self.data['rates']['ioniz']['ionization'][i,:,:],
                                                       self.data['user']['dens_grid']),axis=0)

                self.data['cr_matrix']['cr_loss'][i,i,:,:] = self.data['cr_matrix']['cr'][i,i,:,:] - \
                                      np.sum(np.einsum('ij,k->ijk',self.data['rates']['ioniz']['ionization'][i,:,:],
                                                       self.data['user']['dens_grid']),axis=0)
                
            #level i populating mechanisms
            #these are the transition rates from higher levels into the level i
            self.data['cr_matrix']['cr'][i,0:len(self.data['atomic']['energy']),:,:] = \
                                    self.data['cr_matrix']['cr'][i,0:len(self.data['atomic']['energy']),:,:] + \
                                    self.data['cr_matrix']['A_ji'][:,i,None,None]

            #these are excitation and dexciation into the level i
            self.data['cr_matrix']['cr'][i,0:len(self.data['atomic']['energy']),:,:] = \
                            self.data['cr_matrix']['cr'][i,:len(self.data['atomic']['energy']),:,:] + \
                     np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ij'][:,i,:],
                               self.data['user']['dens_grid']) + \
                     np.einsum('ij,k->ijk',self.data['cr_matrix']['q_ji'][:,i,:],
                               self.data['user']['dens_grid'])


        #recombination out of plus ion to current
        if(self.data['user']['use_recombination']):
            for p in range(0,nsigmaplus):
                #these are the rates from the plus ion into the current ion levels
                self.data['cr_matrix']['cr'][0:len(self.data['atomic']['energy']),
                                             len(self.data['atomic']['energy'])+p,:,:] =\
                self.data['cr_matrix']['cr'][0:len(self.data['atomic']['energy']),
                                             len(self.data['atomic']['energy'])+p,:,:] + \
                np.einsum('ij,k->ijk',self.data['rates']['recomb']['recombination'][:,p,:],
                          self.data['user']['dens_grid'])

                #this is the sum of all the pop lost to the current ion levels taken
                #out of the plus ion
                self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+p,
                                             len(self.data['atomic']['energy'])+p,:,:] = \
                self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+p,
                                             len(self.data['atomic']['energy'])+p,:,:] - \
                np.sum(np.einsum('ij,k->ijk',self.data['rates']['recomb']['recombination'][:,p,:],
                                             self.data['user']['dens_grid']),axis=0)


        if(self.data['user']['use_recombination_three_body']):
            for p in range(0,nsigmaplus):
                self.data['cr_matrix']['cr'][0:len(self.data['atomic']['energy']),
                                             len(self.data['atomic']['energy'])+p,:,:] = \
                self.data['cr_matrix']['cr'][0:len(self.data['atomic']['energy']),
                                             len(self.data['atomic']['energy'])+p,:,:] + \
                np.einsum('ij,k->ijk',self.data['rates']['recomb']['recomb_three_body'][:,p,:],
                                      self.data['user']['dens_grid']**2)
            
                self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+p,
                                             len(self.data['atomic']['energy'])+p,:,:] = \
                self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+p,
                                             len(self.data['atomic']['energy'])+p,:,:] - \
                np.sum(np.einsum('ij,k->ijk',self.data['rates']['recomb']['recomb_three_body'][:,p,:],
                                      self.data['user']['dens_grid']**2),axis=0)
                
        if(self.data['rates']['ioniz']['ionization'].size > 0):                
            for p in range(0,nsigmaplus):
                self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+ p,
                                                 0:len(self.data['atomic']['energy']),:,:] = \
                                     self.data['cr_matrix']['cr'][len(self.data['atomic']['energy'])+p,
                                                         0:len(self.data['atomic']['energy']),:,:] + \
                               np.einsum('ij,k->ijk',self.data['rates']['ioniz']['ionization'][:,p,:],
                                                   self.data['user']['dens_grid'])


        
    def solve_quasi_static(self):
        """This function will solve the CR matrix using the quasistatic approximation after solving this problem
           the generalized radiative coefficients GCRs will be calculated other claculated quanities such as 
           PEC and SXB are also calculated. This function is analgous to ADAS208.
        """
        
        if('cr_matrix' not in self.data.keys()):
            self.populate_cr_matrix()
        levels_to_keep = np.setdiff1d( np.linspace(0,len(self.data['atomic']['energy'])-1,
                                              len(self.data['atomic']['energy']) ,dtype='int64'),
                                       self.data['atomic']['metas'])
        
        #[:,0:len(self.data['atomic']['metas']),:,:] 
        self.data['cr_matrix']['beta']= \
                                          -self.data['cr_matrix']['cr'][levels_to_keep][:,self.data['atomic']['metas']]
        aa_tmp = self.data['cr_matrix']['cr'][levels_to_keep][:,levels_to_keep]

        if(self.data['user']['use_recombination'] or self.data['user']['use_recombination_three_body']):
            #num_recombs = np.shape(self.data['rates']['recomb']['recombination'])[1]
            num_recombs = len(self.data['atomic']['ion_pot'])
            recomb_driving_lvls = len(self.data['atomic']['energy']) + \
                                  np.linspace(0,num_recombs-1,num_recombs,dtype=int)
            self.data['cr_matrix']['beta'] = np.append(self.data['cr_matrix']['beta'],
                                                   -self.data['cr_matrix']['cr'][levels_to_keep][:,recomb_driving_lvls],axis=1)



            

        self.data['cr_matrix']['aa_inv'] = np.zeros((len(aa_tmp),len(aa_tmp), len(self.data['user']['temp_grid']),
                               len(self.data['user']['dens_grid'])))
        '''
        if(self.data['user']['use_ionization']):
            for p in range(0,len(self.data['atomic']['energy'])):
                np.einsum('ij,k->ijk',self.data['rates']['ioniz']['ionization'][p,:,:],
                                                                         self.data['user']['dens_grid'])
        '''
        self.data['cr_matrix']['aa_inv'] = np.linalg.inv(aa_tmp.transpose(2,3,0,1)).transpose(2,3,0,1)

        self.data['processed'] = {}
        self.data['processed']['aa_tmp'] = aa_tmp
        self.data['processed']['excited_levels'] = levels_to_keep
        if(self.data['user']['use_recombination'] or self.data['user']['use_recombination_three_body']):
            self.data['processed']['pops'] = np.zeros((len(aa_tmp),
                                            len(self.data['atomic']['metas'])+len(self.data['atomic']['ion_pot']),
                                            len(self.data['user']['temp_grid']),
                                            len(self.data['user']['dens_grid'])))
        else:
            self.data['processed']['pops'] = np.zeros((len(aa_tmp),
                                            len(self.data['atomic']['metas']),
                                            len(self.data['user']['temp_grid']),
                                            len(self.data['user']['dens_grid'])))
        
        self.data['processed']['pops'] = np.einsum('ijkl,jnkl->inkl',
                                                   self.data['cr_matrix']['aa_inv'],
                                                   self.data['cr_matrix']['beta'])
        #normalization of populations is done here. Quasistatic approximation assumes that
        #there is a small population in 'excited' states, this is sometimes not the case
        #so the populations must be renomalized the default when there is one metastable
        #present is to normalize, but there is a choice for the user to decide
        #default when there is more than one metastable is to not normalize
        if(len(self.data['atomic']['metas']) > 1):

            self.data['processed']['driving_populations_norm'] = False
        else:
            self.data['processed']['driving_populations_norm'] = True

        if(not self.data['user']['default_pop_norm']):
            self.data['processed']['driving_populations_norm'] =\
                          not self.data['processed']['driving_populations_norm']

        if(self.data['processed']['driving_populations_norm']):
            self.data['processed']['pops'] = self.data['processed']['pops'] /\
                                             (1 + np.sum(self.data['processed']['pops'],axis=0))

        #setup pec and wavelength stuff start them as lists because
        #we wont know which levels actually have pec values before we loop through them
        #this lists are eventually converted into np arrays
        self.data['processed']['pecs'] = [] 
        self.data['processed']['pec_levels'] = [] 
        self.data['processed']['wave_vac'] = [] 

        for i in range(0,len(self.data['cr_matrix']['A_ji'])):
            for j in range(0,len(levels_to_keep)):
                #if statement here because there will be a number of transitions that have zero
                #for the PEC value and don't want to have to keep track of these.
                if(levels_to_keep[j] > i and self.data['cr_matrix']['A_ji'][j,i] > 1E-31):
                    self.data['processed']['pecs'].append( self.data['cr_matrix']['A_ji'][levels_to_keep[j],i]*\
                                                        self.data['processed']['pops'][j]/ \
                                                        self.data['user']['dens_grid'])
                    self.data['processed']['wave_vac'].append( 1.e7/abs(self.data['atomic']['energy'][[levels_to_keep[j]]] - \
                                                                                      self.data['atomic']['energy'][i]))
                    self.data['processed']['pec_levels'].append( np.array([levels_to_keep[j],i]))
        
        self.data['processed']['pecs'] = np.asarray(self.data['processed']['pecs'])

        self.data['processed']['wave_vac'] = np.asarray(self.data['processed']['wave_vac'])
        self.data['processed']['wave_air'] = self.data['processed']['wave_vac']/ \
                                                       convert_to_air(self.data['processed']['wave_vac'])
        self.data['processed']['pec_levels'] = np.asarray(self.data['processed']['pec_levels'])

        self.data['processed']['scd'] = np.zeros((len(self.data['atomic']['metas']),
                                                  len(self.data['atomic']['ion_pot']),
                                                  len(self.data['user']['temp_grid']),
                                                  len(self.data['user']['dens_grid'])))
        self.data['processed']['acd'] = np.zeros((len(self.data['atomic']['metas']),
                                                  len(self.data['atomic']['ion_pot']),
                                                  len(self.data['user']['temp_grid']),
                                                  len(self.data['user']['dens_grid'])))
        self.data['processed']['xcd'] = np.zeros((len(self.data['atomic']['ion_pot']),
                                                  len(self.data['atomic']['ion_pot']),
                                                  len(self.data['user']['temp_grid']),
                                                  len(self.data['user']['dens_grid'])))
        self.data['processed']['pop_lvl'] = np.zeros((len(self.data['cr_matrix']['aa_inv']),
                                                       len(self.data['cr_matrix']['aa_inv']),
                                                       len(self.data['atomic']['metas']),
                                                       len(self.data['user']['temp_grid']),
                                                       len(self.data['user']['dens_grid'])))
        #these are how levels get populated
        self.data['processed']['pop_lvl'] = np.einsum('ijkl,jmkl->ijmkl',
                                                      self.data['cr_matrix']['aa_inv'],
                                                      self.data['cr_matrix']['beta'][:,:,:,:])
        #population of levels with no normalization
        self.data['processed']['pops_no_norm'] = np.sum(self.data['processed']['pop_lvl'],axis=1)
        
        if(self.data['processed']['driving_populations_norm']):
            self.data['processed']['pop_lvl'] = self.data['processed']['pop_lvl']/   \
                                            (1+np.sum(np.sum(self.data['processed']['pop_lvl'],axis=1),axis=0))
        #the F matrix
        self.data['processed']['F'] = np.sum(self.data['processed']['pop_lvl'][:,:,0:len(self.data['atomic']['metas']),:,:],axis=1)
        #effective ionization rate
        if(self.data['rates']['ioniz']['ionization'].size > 0):
            self.data['processed']['scd'] = np.einsum('ipk,imkl->mpkl',
                                                      self.data['rates']['ioniz']['ionization'][levels_to_keep,:,:],
                                                      self.data['processed']['F'])

            if(self.data['processed']['driving_populations_norm']):
                self.data['processed']['scd'] = self.data['processed']['scd'] + \
                                              np.einsum('ipk,ikl->ipkl',
                                              self.data['rates']['ioniz']['ionization'][self.data['atomic']['metas'],:,:],
                                              1/(1+np.sum(self.data['processed']['pops_no_norm'],axis=0)))

            else:
                self.data['processed']['scd'] = self.data['processed']['scd'] + \
                                                self.data['rates']['ioniz']['ionization'][self.data['atomic']['metas'],:,:,None]


        if(self.data['user']['use_recombination'] and self.data['user']['use_recombination_three_body']):
            #this is the total recombination with three body and the rates that in included in the adf04 file
            recomb_coeff = np.einsum('ijk,l->ijkl',
                                     self.data['rates']['recomb']['recomb_three_body'],
                                     self.data['user']['dens_grid']) + \
                                     self.data['rates']['recomb']['recombination'][:,:,:,None]
        elif(self.data['user']['use_recombination']):
            recomb_coeff = self.data['rates']['recomb']['recombination'][:,:,:,None]
        elif(self.data['user']['use_recombination_three_body']):            
            recomb_coeff = np.einsum('ijk,l->ijkl',
                                     self.data['rates']['recomb']['recomb_three_body'],
                                     self.data['user']['dens_grid'])
        if(self.data['user']['use_recombination'] or self.data['user']['use_recombination_three_body']):
            #effective recombination rate
            self.data['processed']['acd'] = -np.einsum('njkl,jmkl->nmkl',

                           self.data['cr_matrix']['cr'][np.c_[self.data['atomic']['metas']], levels_to_keep,:,:],
                           np.einsum('ijkl,jmkl->imkl', self.data['cr_matrix']['aa_inv'][0:len(self.data['atomic']['energy']),
                           0:len(self.data['atomic']['energy']),:,:],
                           recomb_coeff[levels_to_keep,:,:,:])
                           )
            
            self.data['processed']['acd'] = self.data['processed']['acd'] + recomb_coeff[self.data['atomic']['metas'],:,:,:]

        self.data['processed']['qcd'] = np.zeros((len(self.data['atomic']['metas']),
                                                  len(self.data['atomic']['metas']),
                                                  len(self.data['user']['temp_grid']),
                                                  len(self.data['user']['dens_grid'])))
        t = []
        for n in range(0,len(self.data['atomic']['metas'])):
            metas_to_keep = np.setdiff1d( self.data['atomic']['metas'],
                                          self.data['atomic']['metas'][n])
            #calculate the metastable cross coupling coefficent, this the number of atoms that start in one
            #metastable and then end up in a different metastable
            metas_to_keep_ind = np.arange(self.data['atomic']['metas'].shape[0])\
                                     [np.in1d(self.data['atomic']['metas'], metas_to_keep)]
            self.data['processed']['qcd'][metas_to_keep_ind,n,:,:] = np.einsum('nkl,l->nkl',
                                                                           
                                    self.data['cr_matrix']['cr'][self.data['atomic']['metas'][n],metas_to_keep,:,:] +\
                                    np.einsum('mkl,mnkl->nkl', self.data['cr_matrix']['cr'][self.data['atomic']['metas'][n],
                                    levels_to_keep,:,:],
                                    self.data['processed']['F'][:,metas_to_keep_ind,:,:]
                                    ),1/self.data['user']['dens_grid'])


        if(self.data['user']['use_recombination'] or self.data['user']['use_recombination_three_body']):            
            for m in range(0,len(self.data['atomic']['ion_pot'])):

                metasplus_to_keep = np.setdiff1d( np.linspace(0,len(self.data['atomic']['ion_pot'])-1,
                                                              len(self.data['atomic']['ion_pot']),dtype='int'),m)
                metasplus_to_keep_ind = np.arange(self.data['atomic']['ion_pot'].shape[0])\
                                     [np.in1d(self.data['atomic']['ion_pot'], metas_to_keep)]
                #parent cross coupling coefficient, start in the parent atom, recombine get redistrubted then
                #ionize back into the parent but into a different parent metastable
                self.data['processed']['xcd'][metasplus_to_keep_ind,np.array([m]),:,:] =-np.einsum('ik,imkl->mkl',

                        self.data['rates']['ioniz']['ionization'][levels_to_keep,m,:],np.einsum('ijkl,jmkl->imkl',
                        self.data['cr_matrix']['aa_inv'][0:len(self.data['atomic']['energy']),0:len(self.data['atomic']['energy']),:,:],
                        recomb_coeff[len(self.data['atomic']['metas']):len(self.data['atomic']['energy']),metasplus_to_keep,:,:])
                        )

    def solve_time_dependent(self):
        """This function will solve the CR matrix with no assumptions. A matrix expoential is used to solve this problem.
           A source term can be included to mimick erosion of fresh atoms or injection of neutral gas or maybe even LIF

           R. LeVeque, Finite Difference Methods for Ordinary and Par-
           tial Differential Equations: Steady-State and Time-Dependent
           Problems (Classics in Applied Mathematics Classics in Applied
           Mathemat), Society for Industrial and Applied Mathematics,
           Philadelphia, PA, USA, 2007.
        """

        if('processed' not in self.data.keys()):
            self.data['processed'] = {}
        self.data['processed']['td'] = {}

        #check if a source term has been provided, if a source term has been provided solve the problem with s asource
        if(self.data['user']['td_source'].any()):

            self.data['processed']['td']['eigenvals'], \
            self.data['processed']['td']['eigenvectors'] = np.linalg.eig(self.data['cr_matrix']['cr'].transpose(2,3,0,1))

            CC = np.dot(np.linalg.inv(self.data['processed']['td']['eigenvectors']),
                                      self.data['user']['td_source'])
            V0 = np.dot(np.linalg.inv(self.data['processed']['td']['eigenvectors']),
                                      self.data['user']['td_n0'])

            eig_zero_ind = np.where(self.data['processed']['td']['eigenvals'] == 0)            
            eig_non_zero = np.delete(self.data['processed']['td']['eigenvals'] ,eig_zero_ind,axis=2)
            amplitude_non = np.delete(V0,eig_zero_ind,axis=2) + np.delete(CC,eig_zero_ind,axis=2)/eig_non_zero
            amplitude_zer = V0[:,:,eig_zero_ind[2]]
            
            v_non = amplitude_non[:,:,:,None]*np.exp(eig_non_zero[:,:,:,None]*self.data['user']['td_t']) - \
                                            np.delete(CC,eig_zero_ind,axis=2)[:,:,:,None]/eig_non_zero[:,:,:,None]
            
            v_zer = CC[:,:,eig_zero_ind[2]][:,:,:,None]*self.data['user']['td_t'] + amplitude_zer[:,:,:,None]
            v = np.insert(v_non,eig_zero_ind[2],v_zer,axis=2)
            self.data['processed']['td']['td_pop'] = np.einsum('klij,kljt->itkl',
                                                                self.data['processed']['td']['eigenvectors'],v)


        else:
            #no source term was provided solve just the normal matrix
            self.data['processed']['td']['eigenvals'], \
            self.data['processed']['td']['eigenvectors'] = np.linalg.eig(self.data['cr_matrix']['cr'].transpose(2,3,0,1))
            
            v0 = np.dot(np.linalg.inv(self.data['processed']['td']['eigenvectors']),self.data['user']['td_n0'])
            vt = v0[:,:,:,None]*np.exp(self.data['processed']['td']['eigenvals'][:,:,:,None]*self.data['user']['td_t'])
            self.data['processed']['td']['td_pop'] = np.einsum('klij,kljt->itkl',
                                                               self.data['processed']['td']['eigenvectors'], vt)
            self.data['processed']['td']['eigenvals'] = self.data['processed']['td']['eigenvals'].transpose(2,0,1)
            self.data['processed']['td']['eigenvectors'] = self.data['processed']['td']['eigenvectors'].transpose(2,3,0,1)


    def split_pec_multiplet(self):
        """This function will solve take LS resolved PECs and split them statistically among the j levels

        """
        if('processed' not in self.data.keys()):
            self.solve_quasi_static()
            
        self.data['processed']['split'] = {}
        self.data['processed']['split']['j_up'] = []
        self.data['processed']['split']['j_low'] = []
        self.data['processed']['split']['pecs'] = []
        self.data['processed']['split']['relative_inten'] = []

        for i in range(0,len(self.data['processed']['pec_levels'])):
            up = self.data['processed']['pec_levels'][i,0]
            low = self.data['processed']['pec_levels'][i,1]
            
            ju,jl,res = split_multiplet( (self.data['atomic']['S'][up]-1)/2.,
                                         self.data['atomic']['L'][up],
                                         (self.data['atomic']['S'][low]-1)/2.,
                                         self.data['atomic']['L'][low])
            
            self.data['processed']['split']['j_low'].append(jl)
            self.data['processed']['split']['j_up'].append(ju)
            self.data['processed']['split']['relative_inten'].append(res)

            if(res.size>0):
                self.data['processed']['split']['pecs'].append(np.einsum('ijk,l->lijk',
                                                                         self.data['processed']['pecs'][i],
                                                                         res/np.sum(res)))
            else:
                self.data['processed']['split']['pecs'].append(self.data['processed']['pecs'][i])
        
    def solve_cr(self):
        if(self.data['user']['use_ionization']):
            self.make_ioniz_from_reduced_ionizrates()
        if(self.data['user']['suppliment_with_ecip']):
            self.make_ecip()
            self.suppliment_with_ecip()

        if(self.data['user']['use_recombination']):
            self.make_recombination_rates_from_file()
        if(self.data['user']['use_recombination_three_body']):
            self.make_three_body_recombination()
        self.make_electron_excitation_rates()
        self.populate_cr_matrix()
        self.solve_quasi_static()
