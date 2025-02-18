import numpy as np
from colradpy.colradpy_class import *
from colradpy.solve_matrix_exponential import *
from scipy.interpolate import RectBivariateSpline
from colradpy.read_adf11 import *

def interp_rates_adf11(logged_temp,logged_dens,temp,dens,logged_gcr): # y are logged_temp and logged_dens args alongside temp and dens? seems redundant?
    # Dane did some optimization, Curt seemed to cobble this together (Dane optimized just this block of code as a self-contained unit)
    # parallelization should be implemented if more performance is needed
    gcr_arr = np.zeros( (np.shape(logged_gcr)[0],np.shape(logged_gcr)[1],len(temp),len(dens)) )

    logged_temp, logged_dens = np.log10(temp), np.log10(dens) # pre-compute to save time
    
    for i in range(0,np.shape(logged_gcr)[0]):
        for j in range(0,np.shape(logged_gcr)[1]):
            interp_gcr = RectBivariateSpline(logged_temp,
                logged_dens,
                logged_gcr[i,j,:,:],
            )
            gcr_arr[i,j] = interp_gcr(logged_temp, logged_dens) # array size matching works when Dane tests it
    return 10**gcr_arr


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

      :param metas: List of arrays for the metastable levels in a charge state
      :type metas: list

      :param temp_grid: Array of the temperatures in (eV)
      :type temp_grid: float array

      :param dens_grid: Array of the densities in (cm-3)
      :type dens_grid: float array

      :param htemp_grid: Temperature grid of thermal CX hydrogen (eV)
      :type htemp_grid: float array

      :param hdens_grid: Density grid of the thermal CX hydrogen densities in (cm-3)
      :type hdens_grid: float array

      :param soln_times: Times to calculate the solution for the time dependent solutions (s)
      :type soln_times: float array

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

      :param use_recombination: Use recombination in the ionization balance. This will produce the ACD and possible XCD rates
                               Can be array of bools or just a single bool,
                             if just a single bool is supplied then all
                             charge states will have the same value.
      :type use_recombination: bool array or bool


      :param use_cx: Use thermal charge exchange  in the ionization balance. 
                     This will produce the CCD rates the input can be array of bools or just a single bool,
                             if just a single bool is supplied then all charge states will have the same value.
      :type use_cx: bool array or bool


      :param keep_species_data: When True this will keep all of the data associated with the CR solutions for each charge state.
                                The can potentially take up lots of memory, the default is False
      :type keep_species_data: bool


      :param init_abund: The initial fractional abundances at t=0
      :type init_abund: float array

      :param source: Source rate into any charge state
      :type source: float array

      :param scale_file_ioniz: Scale ionization in the file
      :type scale_file_ioniz: bool

      :param ne_tau: n_e*tau values for time dependent calc to be evaluated at can not be defined at the same time as td_t (cm-3 s)
      :type ne_tau: float array

      :param adf11_files: the files provided are adf11 formatted files
      :type adf11_files: bool

    """
 
    def __init__(self,fils,metas,temp_grid, dens_grid, htemp_grid=np.array([]), hdens_grid=np.array([]),
                 soln_times=np.array([]), use_ionization=True,
                 suppliment_with_ecip=True, use_recombination_three_body=True,use_recombination=True,
                 use_cx=False,
                 keep_charge_state_data=False,init_abund = np.array([]), source= np.array([]),
                 temp_dens_pair=False,scale_file_ioniz=False, ne_tau=np.array([-1]),adf11_files=False,
                 hdf5_files=False):


        #this is basic input data that is the same whether adf04 or adf11 inputs are used
        self.data = {}        
        self.data['cr_data'] = {}
        self.data['cr_data']['gcrs'] = {}
        self.data['user'] = {}
        self.data['user']['temp_grid'] = np.asarray(temp_grid) #eV
        self.data['user']['dens_grid'] = np.asarray(dens_grid)#cm-3
        self.data['user']['fils'] = np.asarray(fils)
        self.data['user']['init_abund'] = np.asarray(init_abund)
        self.data['user']['soln_times'] = np.asarray(soln_times)





        



        
########################################################################################################################
        if(adf11_files):#adf11 specific formatting
            self.data['input_file'] = {}            
            for i in range(0,len(fils)):
                self.data['user']['temp_dens_pair'] = False
                


                adf11 = read_adf11(fils[i])



                if(type(use_cx) == bool):
                    self.data['user']['use_cx'] = np.ones(len(adf11['input_file']['metas'])-1,dtype=bool)
                    self.data['user']['use_cx'][:] = False
                

                for j in range(0,len(adf11['input_file']['metas'])-1):
                    if( str(j) not in self.data['cr_data']['gcrs']):
                        self.data['cr_data']['gcrs'][str(j)] = {}

                    if( 'scd' in fils[i]):
                        self.data['cr_data']['gcrs'][str(j)]['scd'] = interp_rates_adf11(adf11['input_file']['temp_grid'],
                                                                                         adf11['input_file']['dens_grid'],
                                                                          temp_grid,dens_grid,adf11['input_file'][str(j)])
                        self.data['input_file']['scd'] = adf11['input_file']
                        
                    if( 'acd' in fils[i]):
                        self.data['cr_data']['gcrs'][str(j)]['acd'] = interp_rates_adf11(adf11['input_file']['temp_grid'],
                                                                                         adf11['input_file']['dens_grid'],
                                                                          temp_grid,dens_grid,adf11['input_file'][str(j)])
                        self.data['input_file']['acd'] = adf11['input_file']
                        
                    if( 'qcd' in fils[i]):
                        self.data['cr_data']['gcrs'][str(j)]['qcd']= interp_rates_adf11(adf11['input_file']['temp_grid'],
                                                                                        adf11['input_file']['dens_grid'],
                                                                         temp_grid,dens_grid,adf11['input_file'][str(j)])
                        self.data['input_file']['qcd'] = adf11['input_file']
                        
                    if( 'xcd' in fils[i]):
                        self.data['cr_data']['gcrs'][str(j)]['xcd']= interp_rates_adf11(adf11['input_file']['temp_grid'],
                                                                                        adf11['input_file']['dens_grid'],
                                                                         temp_grid,dens_grid,adf11['input_file'][str(j)])
                        self.data['input_file']['xcd'] = adf11['input_file']                               
                    if( 'ccd' in fils[i]):        
                        self.data['cr_data']['gcrs'][str(j)]['ccd']= interp_rates_adf11(adf11['input_file']['temp_grid'],
                                                                                        adf11['input_file']['dens_grid'],
                                                                         temp_grid,dens_grid,adf11['input_file'][str(j)])
                        self.data['input_file']['ccd'] = adf11['input_file']                        
                        
                    if('qcd' not in self.data['cr_data']['gcrs'][str(j)]):

                        self.data['cr_data']['gcrs'][str(j)]['qcd'] = np.zeros( (np.shape(self.data['cr_data']['gcrs'][str(j)]['scd'])[0],
                                                                                 np.shape(self.data['cr_data']['gcrs'][str(j)]['scd'])[0],
                                                                                 len(temp_grid),
                                                                                 len(dens_grid)) )

                    if('xcd' not in self.data['cr_data']['gcrs'][str(j)]):

                        self.data['cr_data']['gcrs'][str(j)]['xcd'] = np.zeros( (np.shape(self.data['cr_data']['gcrs'][str(j)]['scd'])[0],
                                                                                 np.shape(self.data['cr_data']['gcrs'][str(j)]['scd'])[0],
                                                                                 len(temp_grid),
                                                                                 len(dens_grid)) )


        elif(hdf5_files):

            import h5py
            import hdfdict

            if(type(use_cx) == bool):
                self.data['user']['use_cx'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_cx'][:] = use_cx

            
            self.data['user']['temp_dens_pair'] = False#probably should remove at some point
            self.data['cr_data']['stage_data'] = {}
            for i,j in enumerate(fils):
                c_hdf5 = hdfdict.load(j,lazy=False)
                self.data['cr_data']['gcrs'][str(i)] = {}
                self.data['cr_data']['gcrs'][str(i)]['scd'] = np.copy(c_hdf5['processed']['scd'])
                self.data['cr_data']['gcrs'][str(i)]['acd'] = np.copy(c_hdf5['processed']['acd'])
                self.data['cr_data']['gcrs'][str(i)]['qcd'] = np.copy(c_hdf5['processed']['qcd'])
                self.data['cr_data']['gcrs'][str(i)]['xcd'] = np.copy(c_hdf5['processed']['xcd'])
                self.data['cr_data']['stage_data'][str(i)] = c_hdf5
                
        else:#adf04 specific formating this is the default
            
            # give the option to the user to choose different
            # ionization and recombination settings for each charge state
            # this also just allows the same thing to be specified for all charge states




            self.data['user']['source'] = np.asarray(source)
            self.data['user']['htemp_grid'] = np.asarray(htemp_grid) #eV
            self.data['user']['hdens_grid'] = np.asarray(hdens_grid) #cm-3

            self.data['user']['temp_dens_pair'] = temp_dens_pair
            self.data['user']['ne_tau'] = ne_tau
            self.data['user']['metas'] = np.asarray(metas)        
            if(self.data['user']['init_abund'].size <1):
                self.data['user']['init_abund'] = np.zeros((np.sum(list(map(len,metas))) +1))
                self.data['user']['init_abund'][0] = 1.


            if( (self.data['user']['ne_tau'][0] != -1) or (self.data['user']['ne_tau'][0] != -1)):
                if(self.data['user']['ne_tau'][0] == -1):
                    nt = ne_tau
                else:
                    nt = self.data['user']['ne_tau']
                td_t = np.einsum('i,j->ij', nt, 1/self.data['user']['dens_grid'])
                print('ne_tau was not equal to -1 so overwiting solution times with ne_tau calculated values')
                self.data['user']['soln_times'] = td_t 





            
            if(type(use_ionization) == bool):
                self.data['user']['use_ionization'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_ionization'][:] = use_ionization
            else:
                self.data['user']['use_ionization'] = use_ionization
            if(type(suppliment_with_ecip) == bool):
                self.data['user']['suppliment_with_ecip'] = np.ones_like(fils,dtype=bool)
                self.data['user']['suppliment_with_ecip'][:] = suppliment_with_ecip
            else:
                self.data['user']['suppliment_with_ecip'] = suppliment_with_ecip
            if(type(use_recombination_three_body) == bool):
                self.data['user']['use_recombination_three_body'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_recombination_three_body'][:] = use_recombination_three_body
            else:
                self.data['user']['use_recombination_three_body'] = use_recombination_three_body
            if(type(use_recombination) == bool):
                self.data['user']['use_recombination'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_recombination'][:] = use_recombination
            else:
                self.data['user']['use_recombination'] = use_recombination
            if(type(use_cx) == bool):
                self.data['user']['use_cx'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_cx'][:] = use_cx
            else:
                self.data['user']['use_cx'] = use_cx

            if(type(use_ionization) == bool):
                self.data['user']['use_ionization'] = np.ones_like(fils,dtype=bool)
                self.data['user']['use_ionization'][:] = use_ionization
            else:
                self.data['user']['use_ionization'] = use_ionization
            if(type(scale_file_ioniz) == bool):
                self.data['user']['scale_file_ioniz'] = np.ones_like(fils,dtype=bool)
                self.data['user']['scale_file_ioniz'][:] = scale_file_ioniz
            else:
                self.data['user']['scale_file_ioniz'] = scale_file_ioniz

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
                print(i)
                tmp = colradpy(fil=str(fils[i]),metas=meta,temp_grid=temp_grid,electron_den=dens_grid,
                               htemp_grid = htemp_grid, hdens_grid = hdens_grid,
                               use_ionization = self.data['user']['use_ionization'][i],
                               suppliment_with_ecip = self.data['user']['suppliment_with_ecip'][i],
                               use_recombination_three_body = self.data['user']['use_recombination_three_body'][i],
                               use_recombination = self.data['user']['use_recombination'][i],
                               use_cx = self.data['user']['use_cx'][i],
                               temp_dens_pair = self.data['user']['temp_dens_pair'])

                tmp.solve_cr() #Solving the CR set of equations
                #keep all the CR data if requested
                if(keep_charge_state_data):
                    self.data['cr_data']['stage_data'][str(i)] = tmp.data
                #keeping the GCR data, ionization can be run from this data only
                self.data['cr_data']['gcrs'][str(i)] = {}
                self.data['cr_data']['gcrs'][str(i)]['qcd'] = tmp.data['processed']['qcd']
                self.data['cr_data']['gcrs'][str(i)]['scd'] = tmp.data['processed']['scd']
                self.data['cr_data']['gcrs'][str(i)]['acd'] = tmp.data['processed']['acd']
                self.data['cr_data']['gcrs'][str(i)]['xcd'] = tmp.data['processed']['xcd']
                if(self.data['user']['use_cx'][i]):
                    self.data['cr_data']['gcrs'][str(i)]['ccd'] = tmp.data['processed']['ccd']

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



        if(self.data['user']['temp_dens_pair']):
            self.data['ion_matrix'] = np.zeros((m, m,len(self.data['user']['temp_grid'])))


            m = 0
            for i in range(0,len(self.data['cr_data']['gcrs'])):
                num_met = np.shape(self.data['cr_data']['gcrs'][str(i)]['qcd'])[0]
                num_ion = np.shape(self.data['cr_data']['gcrs'][str(i)]['scd'])[1]
                diag_met = np.diag_indices(num_met)
                diag_ion = np.diag_indices(num_ion)

                #populate QCDs in ion balance
                self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] \
                                       -np.sum(self.data['cr_data']['gcrs'][str(i)]['qcd'],axis=1)

                self.data['ion_matrix'][m:m+num_met,m:m+num_met,:]=self.data['ion_matrix'][m:m+num_met,m:m+num_met,:]\
                                        +self.data['cr_data']['gcrs'][str(i)]['qcd'].transpose(1,0,2)

                #populate SCDs in ion balance
                self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]] = self.data['ion_matrix'][m+diag_met[0],m+diag_met[1]]\
                                          -np.sum(self.data['cr_data']['gcrs'][str(i)]['scd'],axis=1)

                self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:] = \
                                    self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m:m+num_met,:] \
                                                       +self.data['cr_data']['gcrs'][str(i)]['scd'].transpose(1,0,2)

                #populate ACDs in ion balance
                self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                                   np.sum(self.data['cr_data']['gcrs'][str(i)]['acd'],axis=0)

                self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:] = \
                                self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:] \
                                                 +  self.data['cr_data']['gcrs'][str(i)]['acd']

                #populate CCDs in ion balance



                
                if(self.data['user']['use_cx'][i]):
                    self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                                self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                                       np.sum(np.einsum('ijn,n->ijn',self.data['cr_data']['gcrs'][str(i)]['ccd'],
                                                                        self.data['user']['hdens_grid']/self.data['user']['dens_grid']),axis=0)


                    self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:] = \
                                    self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:] \
                                                     +  np.einsum('ijn,n->ijn',self.data['cr_data']['gcrs'][str(i)]['ccd'],
                                                                  self.data['user']['hdens_grid']/self.data['user']['dens_grid'])                                                                  

                #populate XCD in ion balance
                self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                            np.sum(self.data['cr_data']['gcrs'][str(i)]['xcd'],axis=1)

                self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:]=\
                            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:] + \
                            self.data['cr_data']['gcrs'][str(i)]['xcd'].transpose(1,0,2)

                m = m + num_met

            
        else:
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
                #populate CCDs in ion balance
                if(self.data['user']['use_cx'][i]):                
                    self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                                   np.sum(np.einsum('ijkn,n->ijkn',self.data['cr_data']['gcrs'][str(i)]['ccd'],
                                                                    self.data['user']['hdens_grid']/self.data['user']['dens_grid']),axis=0)

                    self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] = \
                                    self.data['ion_matrix'][m:m+num_met,m+num_met:m+num_met+num_ion,:,:] \
                                                     +  np.einsum('ijkn,n->ijkn',self.data['cr_data']['gcrs'][str(i)]['ccd'],
                                                                  self.data['user']['hdens_grid']/self.data['user']['dens_grid'])

                #populate XCD in ion balance
                self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                            self.data['ion_matrix'][m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                            np.sum(self.data['cr_data']['gcrs'][str(i)]['xcd'],axis=1)

                self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:]=\
                            self.data['ion_matrix'][m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:] + \
                            self.data['cr_data']['gcrs'][str(i)]['xcd'].transpose(1,0,2,3)

                m = m + num_met


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
            self.data['processed']['pops_td'] = 0
            self.data['processed']['eigen_val'] = 0
            self.data['processed']['eigen_vec'] = 0            
            #solve the ionization balance set of equation with no source term
        if(self.data['user']['temp_dens_pair']):
            self.data['processed']['pops_td'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential(
                                               np.einsum('ijk,k->ijk',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']),n0,td_t)

        else:
            self.data['processed']['pops_td'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential(
                                               np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']),n0,td_t)
    def solve_source(self,n0=np.array([]), s0=np.array([]), td_t=np.array([]),ne_tau = np.array([-1])):
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

        if( (self.data['user']['ne_tau'][0] != -1) or (self.data['user']['ne_tau'][0] != -1)):
            if(self.data['user']['ne_tau'][0] == -1):
                nt = ne_tau
            else:
                nt = self.data['user']['ne_tau']
            td_t = np.einsum('i,j->ij', nt, 1/self.data['user']['dens_grid'])
            print('ne_tau was not equal to -1 so overwiting solution times with ne_tau calculated values')
            self.data['user']['soln_times'] = td_t 
            
        if(s0.size < 1):
            if(self.data['user']['source'].size < 1):
                s0 = np.zeros((np.sum(list(map(len,self.data['user']['metas']))) +1))
            else:
                s0 = self.data['user']['source']


                
        if('processed' not in self.data.keys()):
            self.data['processed'] = {}




        if(self.data['user']['temp_dens_pair']):
            self.data['processed']['pops_td_source'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential_source(
                                               np.einsum('ijk,k->ijk',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']),n0,s0,td_t)
        else:
            self.data['processed']['pops_td_source'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential_source(
                                               np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']),n0,s0,td_t)

        self.data['processed']['pops_source'] = self.data['processed']['pops_td_source']

        '''
        self.data['processed']['pops_td'],\
        self.data['processed']['eigen_val'],\
        self.data['processed']['eigen_vec'] = solve_matrix_exponential_source(
                                                 np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                                 self.data['user']['dens_grid']),n0,s0,td_t)
        '''

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
        if(self.data['user']['temp_dens_pair']):
            self.data['processed']['pops_ss'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential_steady_state(
                                               np.einsum('ijk,k->ijk',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']))
        else:
            self.data['processed']['pops_ss'],\
                self.data['processed']['eigen_val'],\
                self.data['processed']['eigen_vec'] = solve_matrix_exponential_steady_state(
                                               np.einsum('ijkl,l->ijkl',self.data['ion_matrix'],
                                               self.data['user']['dens_grid']))
