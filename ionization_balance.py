import numpy as np
from colradpy_class import colradpy
from solve_matrix_exponential import *

class ionization_balance():

    def __init__(self, fils,metas,temp_grid, dens_grid,use_ionization=True,
                suppliment_with_ecip=True, use_recombination_three_body=True,use_recombination=True):

        self.ion = []
        self.metas_tot = 0
        self.temp_grid = np.asarray(temp_grid)
        self.dens_grid = np.asarray(dens_grid)
        for i in range(0,len(fils)):
            self.ion.append(colradpy(fils[i],metas[i],temp_grid,dens_grid,use_ionization=use_ionization,
                                     suppliment_with_ecip=suppliment_with_ecip,
                                     use_recombination_three_body=use_recombination_three_body,
                                     use_recombination=use_recombination))
            self.ion[i].solve_cr()
            #self.metas_tot = self.metas_tot + self.ion[i].data['atomic']['metas']
        self.populate_ion_matrix()

    def populate_ion_matrix(self):

        m = 0
        for i in range(0,len(self.ion)):
            m = m + len(self.ion[i].data['atomic']['metas'])
            if(i == len(self.ion) -1):
                m = m+len(self.ion[i].data['atomic']['ion_pot'])
        self.ion_matrix = np.zeros((m, m,len(self.temp_grid), len(self.dens_grid)))

        m = 0
        for i in range(0,len(self.ion)):
            num_met = len(self.ion[i].data['atomic']['metas'])
            num_ion = len(self.ion[i].data['atomic']['ion_pot'])
            diag_met = np.diag_indices(len(self.ion[i].data['atomic']['metas']))
            diag_ion = np.diag_indices(num_ion)

            #populate QCDs in ion balance
            self.ion_matrix[m+diag_met[0],m+diag_met[1]] = self.ion_matrix[m+diag_met[0],m+diag_met[1]] \
                                   -np.sum(self.ion[i].data['processed']['qcd'],axis=1)
            
            self.ion_matrix[m:m+num_met,m:m+num_met,:,:]=self.ion_matrix[m:m+num_met,m:m+num_met,:,:]\
                                    +self.ion[i].data['processed']['qcd'].transpose(1,0,2,3)

            #populate SCDs in ion balance

            self.ion_matrix[m+diag_met[0],m+diag_met[1]] = self.ion_matrix[m+diag_met[0],m+diag_met[1]]\
                                      -np.sum(self.ion[i].data['processed']['scd'],axis=1)
            
            self.ion_matrix[m+num_met:m+num_met+num_ion,m:m+num_met,:,:] = \
                                self.ion_matrix[m+num_met:m+num_met+num_ion,m:m+num_met,:,:] \
                                                   +self.ion[i].data['processed']['scd'].transpose(1,0,2,3)

            #populate ACDs in ion balance
            self.ion_matrix[m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.ion_matrix[m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                               np.sum(self.ion[i].data['processed']['acd'],axis=0)

            self.ion_matrix[m:m+num_met,m+num_met:m+num_met+num_ion,:,:] = \
                            self.ion_matrix[m:m+num_met,m+num_met:m+num_met+num_ion,:,:] \
                                             +  self.ion[i].data['processed']['acd']

            #populate XCD in ion balance
            self.ion_matrix[m+num_met+diag_ion[0], m+num_met+diag_ion[1] ] = \
                                        self.ion_matrix[m+num_met+diag_ion[0], m+num_met+diag_ion[1] ]-\
                                        np.sum(self.ion[i].data['processed']['xcd'],axis=1)

            self.ion_matrix[m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:]=\
                        self.ion_matrix[m+num_met:m+num_met+num_ion,m+num_met:m+num_met+num_ion,:,:] + \
                        self.ion[i].data['processed']['xcd'].transpose(1,0,2,3)
            
            m = m + num_met
            #self.ion_matrix[m+num_met:m+num_met+num_ion,m:m+num_met,:,:]\

    def solve_no_source(self,n0,td_t):
        self.pops,self.eigen_val,self.eigen_vec = solve_matrix_exponential(
                                       np.einsum('ijkl,l->ijkl',self.ion_matrix,self.dens_grid),n0,td_t)

    def solve_source(self,n0,s0,td_t):
        self.pops,self.eigen_val,self.eigen_vec = solve_matrix_exponential_source(
                                       np.einsum('ijkl,l->ijkl',self.ion_matrix,self.dens_grid),n0,s0,td_t)

