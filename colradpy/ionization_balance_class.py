import numpy as np
from .colradpy_class import colradpy
from .solve_matrix_exponential import *

class ionization_balance():
    """The ionization balance class takes adf04 file inputs runs CR calcualations
       'colradpy_class.py' to get GCR coefficients. The GCR coefficients are then
       assemble into the ionization balance matrix. This matrix is solved in the 
       same way as the time dependent CR problem with 'solve_matrix_exponential.py'.
       The matrix can be solved time dependently with and without a source.
       The time-independent solution can also be solved. See dictionary structure documentation.

    :Args:       
      :param fils: Array of the input adf04 file
      :type fils: string array

      :param temp_grid: Array of the temperatures in eV
      :type temp_grid: float array

      :param dens_grid: Array of the densities in cm-3
      :type dens_grid: float array

      :param use_ionization: Choose if ionization will be use. This will proced the SCD rates.
       This should probably always be true if running an ionization balance but the option is there.
                             Can be array of bools or just a single bool,
                             if just a single bool is supplied then all
                             charge states will have the same value.

      :type use_ionization: bool array or bool

      :param suppliment_with_ecip: Suppliments ionization rates with ECIP when ionization rates are not 
                                   present in the adf04 file.
                               Can be array of bools or just a single bool,
                             if just a single bool is supplied then all
                             charge states will have the same value.

      :type suppliment_with_ecip: bool array or bool

      :param use_recombination_three_body: Adds in three body recombination, there must be ionization present for this to work
                               Can be array of bools or just a single bool,
                             if just a single bool is supplied then all
                             charge states will have the same value.

      :type use_recombination_three_body: bool array or bool

      :param use_recombination: Use recombination in the ionization balance. THis will produce the ACD and possible XCD rates
                               Can be array of bools or just a single bool,
                             if just a single bool is supplied then all
                             charge states will have the same value.

      :type use_recombination: bool array or bool


      :param keep_species_data: When True this will keep all of the data associated with the CR solutions for each charge state.
                                The can potentially take up lots of memory, the default is False

      :type use_recombination: bool


    """
 
    def __init__(self,fils,metas,temp_grid, dens_grid, soln_times=np.array([]), use_ionization=True,
                 suppliment_with_ecip=True, use_recombination_three_body=True,use_recombination=True,
                 keep_charge_state_data=False,init_abund = np.array([]), source= np.array([])):
        self.data = {}

        self.data['cr_data'] = {}
        self.data['cr_data']['gcrs'] = {}
        self.data['user'] = {}
        self.data['user']['temp_grid'] = np.asarray(temp_grid) #eV
        self.data['user']['dens_grid'] = np.asarray(dens_grid)#cm-3
        self.data['user']['fils'] = np.asarray(fils)
        self.data['user']['init_abund'] = np.asarray(init_abund)
        self.data['user']['soln_times'] = np.asarray(soln_times)
        self.data['user']['source'] = np.asarray(source)
        self.data['user']['metas'] = np.asarray(metas)
        
        if(self.data['user']['init_abund'].size <1):
            self.data['user']['init_abund'] = np.zeros((np.sum(list(map(len,metas))) +1))
            self.data['user']['init_abund'][0] = 1.
        
        # give the option to the user to choose different
        # ionization and recombination settings for each charge state
        # this also just allows the same thing to be specified for all charge states
        if(type(use_ionization) == bool):
            self.data['user']['use_ionization'] = np.ones_like(fils,dtype=bool)
            self.data['user']['use_ionization'][:] = use_ionization

        if(type(suppliment_with_ecip) == bool):
            self.data['user']['suppliment_with_ecip'] = np.ones_like(fils,dtype=bool)
            self.data['user']['suppliment_with_ecip'][:] = suppliment_with_ecip

        if(type(use_recombination_three_body) == bool):
            self.data['user']['use_recombination_three_body'] = np.ones_like(fils,dtype=bool)
            self.data['user']['use_recombination_three_body'][:] = use_recombination_three_body

        if(type(use_recombination) == bool):
            self.data['user']['use_recombination'] = np.ones_like(fils,dtype=bool)
            self.data['user']['use_recombination'][:] = use_recombination

        if(type(use_ionization) == bool):
            self.data['user']['use_ionization'] = np.ones_like(fils,dtype=bool)
            self.data['user']['use_ionization'][:] = use_ionization


        self.data['user']['keep_charge_state_data'] = keep_charge_state_data

        #default is to not keep all of the data from the CR calculations
        #this can take up a lot of memory and the general use case for this class only
        #calls for the use of GCR coefficients so its kind of a waste to keep all the data
        if(keep_charge_state_data):
            self.data['cr_data']['stage_data'] = {}
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

            tmp = colradpy(str(fils[i]),meta,temp_grid,dens_grid,
                           self.data['user']['use_ionization'][i],
                           self.data['user']['suppliment_with_ecip'][i],
                           self.data['user']['use_recombination_three_body'][i],
                           self.data['user']['use_recombination'][i])
            
            tmp.solve_cr()
            #keep all the CR data if requested
            if(keep_charge_state_data):
                self.data['cr_data']['stage_data'][str(i)] = tmp.data
            #keeping the GCR data, ionization can be run from this data only
            self.data['cr_data']['gcrs'][str(i)] = {}
            self.data['cr_data']['gcrs'][str(i)]['qcd'] = tmp.data['processed']['qcd']
            self.data['cr_data']['gcrs'][str(i)]['scd'] = tmp.data['processed']['scd']
            self.data['cr_data']['gcrs'][str(i)]['acd'] = tmp.data['processed']['acd']
            self.data['cr_data']['gcrs'][str(i)]['xcd'] = tmp.data['processed']['xcd']            

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
        m = np.shape(self.data['cr_data']['gcrs']['0']['qcd'])[0]#num metastables in ground
        m = m + np.shape(self.data['cr_data']['gcrs']['0']['scd'])[1]#num metastables in plus

        for i in range(1,len(self.data['cr_data']['gcrs'])):
            m = m + np.shape(self.data['cr_data']['gcrs'][str(i)]['scd'])[1]#num of metastales in the plus loop
        #create an empty matrix to hold the GCR rates
        self.data['ion_matrix'] = np.zeros((m, m,len(self.data['user']['temp_grid']),
                                                 len(self.data['user']['dens_grid'])))

        m = 0
        for i in range(0,len(self.data['cr_data']['gcrs'])):
            num_met = np.shape(self.data['cr_data']['gcrs'][str(i)]['qcd'])[0]
            num_ion = np.shape(self.data['cr_data']['gcrs'][str(i)]['scd'])[1]
            diag_met = np.diag_indices(num_met)
            diag_ion = np.diag_indices(num_ion)

            #populate QCDs in ion balance
            self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] \
                                   -np.sum(self.data['cr_data']['gcrs'][str(i)]['qcd'],axis=1)

            self.data['ion_matrix'][m:m+num_met,m:m+num_met,:,:]=self.data['ion_matrix'][m:m+num_met,m:m+num_met,:,:]\
                                    +self.data['cr_data']['gcrs'][str(i)]['qcd'].transpose(1,0,2,3)

            #populate SCDs in ion balance
            self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]]\
                                      -np.sum(self.data['cr_data']['gcrs'][str(i)]['scd'],axis=1)

            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:] = \
                                self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:] \
                                                   +self.data['cr_data']['gcrs'][str(i)]['scd'].transpose(1,0,2,3)

            #populate ACDs in ion balance
            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                               np.sum(self.data['cr_data']['gcrs'][str(i)]['acd'],axis=0)

            self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] = \
                            self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] \
                                             +  self.data['cr_data']['gcrs'][str(i)]['acd']

            #populate XCD in ion balance
            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                        np.sum(self.data['cr_data']['gcrs'][str(i)]['xcd'],axis=1)

            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:]=\
                        self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:] + \
                        self.data['cr_data']['gcrs'][str(i)]['xcd'].transpose(1,0,2,3)

            m = m + num_met
            #self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:,:]\

    def solve_no_source(self,n0=np.array([]),td_t=np.array([])):
        """Solves the ionization balance matrix given an initial populations and times

        :param n0: The initial populations of levels at t=0
        :type n0: float array

        :param td_t: Solution times
        :type td_t: float array


        populates the pops, eigen_val and eigen_vec

        """

        if(n0.size < 1):
            n0 = self.data['user']['init_abund']
        if(td_t.size < 1):
            td_t = self.data['user']['soln_times']

        if('processed' not in self.data.keys()):
            self.data['processed'] = {}
            self.data['processed']['pops'] = 0
            self.data['processed']['eigen_val'] = 0
            self.data['processed']['eigen_vec'] = 0

            self.data['processed']['pops'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential(
                                               np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']),n0,td_t)

    def solve_source(self,n0=np.array([]), s0=np.array([]), td_t=np.array([])):
        """Solves the ionization balance matrix given an
             initial populations and times when a source term is present

        :param n0: The initial populations of levels at t=0
        :type n0: float array

        :param td_t: Solution times
        :type td_t: float array


        populates the pops, eigen_val and eigen_vec

        """

        if(n0.size < 1):
            n0 = self.data['user']['init_abund']
        if(td_t.size < 1):
            td_t = self.data['user']['soln_times']
        if(s0.size < 1):
            if(self.data['user']['source'].size < 1):
                s0 = np.zeros((np.sum(list(map(len,self.data['user']['metas']))) +1))
            else:
                s0 = self.data['user']['source']



                
        if('processed' not in self.data.keys()):
            self.data['processed'] = {}

        self.data['processed']['pops'],\
        self.data['processed']['eigen_val'],\
        self.data['processed']['eigen_vec'] = solve_matrix_exponential_source(
                                                 np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                                 self.data['user']['dens_grid']),n0,s0,td_t)



    def solve_time_independent(self):
        """Solves the ionization balance matrix for the steady-state (time-independent) solution.

           This is going to use the time dependent method just solving at 8 e-folding times
           for the second to smallest eigen value. Note that there are faster methods to do this
           but its more work to code it up and the time-dependent part is already done.
           This function is mostly a niceity for the casual user.
           The smallest eigenvalue corresponds to the
           steady state so we want to wait until the second to last componet completely dies off.
           This is done for the smallest over the given temperature and density parameter space
           this has been tested for 17 orders of magnitude and I haven't run into overflow problems.


        :param n0: The initial populations of levels at t=0
        :type n0: float array

        :param td_t: Solution times
        :type td_t: float array

        populates the pops, eigen_val and eigen_vec
        """
        if('processed' not in self.data.keys()):
            self.data['processed'] = {}
            self.data['processed']['pops_ss'] = 0
            
        #for the steady state the initial conditions don't matter so just make it simple        
        n0 = np.zeros((len(self.data['ion_matrix'])))
        n0[0] = 1.

        #run the first time at t=0 just to get eigenvals
        self.data['processed']['pops_ss'],\
        self.data['processed']['eigen_val'],\
        self.data['processed']['eigen_vec'] = solve_matrix_exponential(
                                                np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                                self.data['user']['dens_grid']),n0,np.array([0]))
        #get the second smallest time
        time = np.max(8/np.abs(np.sort(self.data['processed']['eigen_val'][:,:,:],
                                       axis=0)[len(self.data['ion_matrix']) -2]))

        #solve at the 8 e-fold time
        self.data['processed']['pops_ss'],\
        self.data['processed']['eigen_val'],\
        self.data['processed']['eigen_vec'] = solve_matrix_exponential(
                                                np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                                self.data['user']['dens_grid']),n0,np.array([time]))
