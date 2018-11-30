import numpy as np
from colradpy_class import colradpy

class ionization_balance():

    def __init__(self, fils,temp_grid, dens_grid):

        self.ion = []
        self.metas_tot = 0
        self.temp_grid = temp_grid
        self.dens_grid = dens_grid
        for i in range(0,len(fils)):
            self.ion.append(colradpy(fils[i],[0,1,2],temp_grid,dens_grid,use_ionization=True,
                suppliment_with_ecip=True, use_recombination_three_body=True,use_recombination=True))
            self.ion[i].solve_cr()
            self.metas_tot = self.metas_tot + self.ion[i].data['atomic']['metas']
        self.populate_ion_matrix()
    def populate_ion_matrix(self):
        
        self.ion_matrix = np.zeros((4, 4,
                                    len(self.temp_grid), len(self.dens_grid)))

        m = 0
        for i in range(0,len(self.ion)):
            num_met = len(self.ion[i].data['atomic']['metas'])
            num_ion = len(self.ion[i].data['atomic']['ion_pot'])
            diag_met = np.diag_indices(len(self.ion[i].data['atomic']['metas']))
            diag_ion = np.diag_indices(num_ion)


            #populate QCDs in ion balance
            self.ion_matrix[diag_met] = self.ion_matrix[diag_met] \
                                   -np.sum(self.ion[i].data['processed']['qcd'],axis=1)
            self.ion_matrix[m:m+num_met,m:m+num_met,:,:]=self.ion_matrix[m:m+num_met,m:m+num_met,:,:]\
                                    +self.ion[i].data['processed']['qcd']
            #populate SCDs in ion balance
            self.ion_matrix[diag_met] = self.ion_matrix[diag_met]\
                                      -np.sum(self.ion[i].data['processed']['scd'],axis=1)
            self.ion_matrix[m:m+num_met,m+num_met:m+num_met+num_ion,:,:] = \
                                         self.ion_matrix[m:m+num_met,m+num_met:m+num_met+num_ion,:,:]\
                                                   +self.ion[i].data['processed']['scd']

            #populate ACDs in ion balance


            
            self.ion_matrix[m+num_met:m+num_met+num_ion,m:m+num_met,:,:] = 
            
