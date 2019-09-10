import numpy as np
from .colradpy_class import colradpy
from .solve_matrix_exponential import *

class ionization_balance():

    def __init__(self,fils,metas,temp_grid, dens_grid,use_ionization=True,
                 suppliment_with_ecip=True, use_recombination_three_body=True,use_recombination=True,
                 keep_species_data=False):

        """ Sets up the dictionaries to hold the data.
            
        """
        self.data = {}
        self.data['gcrs'] = {}
        
        self.data['user'] = {}
        self.data['user']['temp_grid'] = np.asarray(temp_grid) #eV
        self.data['user']['dens_grid'] = np.asarray(dens_grid)#cm-3
        self.data['user']['use_ionization'] = use_ionization
        self.data['user']['suppliment_with_ecip'] = suppliment_with_ecip
        self.data['user']['use_recombination_three_body'] = use_recombination_three_body
        self.data['user']['use_recombination'] = use_recombination
        self.data['user']['keep_species_data'] = keep_species_data
        
        #default is to not keep all of the data from the CR calculations
        #this can take up a lot of memory and the general use case for this class only
        #calls for the use of GCR coefficients so its kind of a waste to keep all the data
        if(keep_species_data):
            self.data['stage_data'] = {}
        #cycle over all of the ion stages that the user chose and run the CR calculation
        #to get the GCR coefficients for that stage
        for i in range(0,len(fils)):
            #default is for user to just give metastables of the ground state and then
            #choose metastables just based off of ionization potentials in adf04 file
            #the user can also give a list of metastable states for every ionstage
            #and this will override the metastables in the adf04 files
            if(type(metas) == list):
                meta = metas[i]
            else:
                if(i == 0):
                    meta = metas
                else:
                    m = np.shape(self.data['gcrs'][str(i -1)]['scd'])[1] 
                    meta = np.linspace(0,m-1,m,dtype=int)
            #setup the CR

            tmp = colradpy(str(fils[i]),meta,temp_grid,dens_grid,use_ionization,
                              suppliment_with_ecip,use_recombination_three_body,use_recombination)
            tmp.solve_cr()
            #keep the data if requested
            if(keep_species_data):
                self.data['stage_data'][str(i)] = tmp.data
            #keeping all of the GCR data, ionization can be run from this data only
            self.data['gcrs'][str(i)] = {}
            self.data['gcrs'][str(i)]['qcd'] = tmp.data['processed']['qcd']
            self.data['gcrs'][str(i)]['scd'] = tmp.data['processed']['scd']
            self.data['gcrs'][str(i)]['acd'] = tmp.data['processed']['acd']
            self.data['gcrs'][str(i)]['xcd'] = tmp.data['processed']['xcd']            

    def populate_ion_matrix(self):
        """This will populate the ionization matrix from the various GCR coefficients
           for all of the ionization stages
           
           This matrix in mostly zeros because we don't have ways to connect charge states
           that are more than one change away. For example if only ground states are included
           this matrix becomes diagonal.

           QCDs bring the atom between metastables states in a charge state
           SCDS bring the atom to the next charge state
           ACDS bring the atom to the previous charge state
           XCDS bring the atom to next charge state between the metastables of that level,
           through the previous charge state

           The columns of this matrix sum to zero
        """
        #finding the total number of states to be tracked (metastables)
        m = np.shape(self.data['gcrs']['0']['qcd'])[0]#num metastables in ground
        m = m + np.shape(self.data['gcrs']['0']['scd'])[1]#num metastables in plus
        
        for i in range(1,len(self.data['gcrs'])):
            m = m + np.shape(self.data['gcrs'][str(i)]['scd'])[1]#num of metastales in the plus loop
        #create an empty matrix to hold the GCR rates
        self.data['ion_matrix'] = np.zeros((m, m,len(self.data['user']['temp_grid']), len(self.data['user']['dens_grid'])))

        m = 0
        for i in range(0,len(self.data['gcrs'])):
            num_met = np.shape(self.data['gcrs'][str(i)]['qcd'])[0]
            num_ion = np.shape(self.data['gcrs'][str(i)]['scd'])[1]
            diag_met = np.diag_indices(num_met)
            diag_ion = np.diag_indices(num_ion)

            #populate QCDs in ion balance
            self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] \
                                   -np.sum(self.data['gcrs'][str(i)]['qcd'],axis=1)
            
            self.data['ion_matrix'][m:m+num_met,m:m+num_met,:,:]=self.data['ion_matrix'][m:m+num_met,m:m+num_met,:,:]\
                                    +self.data['gcrs'][str(i)]['qcd'].transpose(1,0,2,3)

            #populate SCDs in ion balance
            self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]]\
                                      -np.sum(self.data['gcrs'][str(i)]['scd'],axis=1)
            
            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:] = \
                                self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:] \
                                                   +self.data['gcrs'][str(i)]['scd'].transpose(1,0,2,3)

            #populate ACDs in ion balance
            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                               np.sum(self.data['gcrs'][str(i)]['acd'],axis=0)

            self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] = \
                            self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] \
                                             +  self.data['gcrs'][str(i)]['acd']

            #populate XCD in ion balance
            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                        np.sum(self.data['gcrs'][str(i)]['xcd'],axis=1)

            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:]=\
                        self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:] + \
                        self.data['gcrs'][str(i)]['xcd'].transpose(1,0,2,3)
            
            m = m + num_met
            #self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:]\

    def solve_no_source(self,n0,td_t):
        """Solves the ionization balance matrix given an initial populations and times

        :param n0: The initial populations of levels at t=0
        :type n0: float array

        :param td_t: Solution times
        :type td_t: float array


        populates the pops, eigen_val and eigen_vec
        
        """
        self.data['pops'],self.data['eigen_val'],self.data['eigen_vec'] = solve_matrix_exponential(
                                       np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],self.data['user']['dens_grid']),n0,td_t)

    def solve_source(self,n0,s0,td_t):
        """Solves the ionization balance matrix given an initial populations and times when a source term is present

        :param n0: The initial populations of levels at t=0
        :type n0: float array

        :param td_t: Solution times
        :type td_t: float array


        populates the pops, eigen_val and eigen_vec
        
        """
        self.data['pops'],self.data['eigen_val'],self.data['eigen_vec'] = solve_matrix_exponential_source(
                                       np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],self.data['user']['dens_grid']),n0,s0,td_t)

