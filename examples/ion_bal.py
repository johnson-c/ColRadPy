import sys
sys.path.append('../')
from colradpy_class import colradpy
import numpy as np
import matplotlib.pyplot as plt
from ionization_balance_class import ionization_balance

#the adf04 files
fils = np.array(['cpb03_ls#be0.dat','cpb03_ls#be1.dat','be2_adf04','be3_adf04'])
temp = np.linspace(5,100,5) #temp grid
dens = np.array([1.e11,1.e14]) #density grid
metas = [np.array([0,1]),np.array([0]),np.array([0,1]),np.array([0])]#number of metastable
                           #this should match the metastables at the top of the adf04 file
                           #this information is used to calculate the QCD values
                           #without it only the SCD, ACD and XCD for a species will be calculated

time = np.linspace(0,.01,1.e4)
                           
ion = ionization_balance(fils, metas, temp, dens, keep_species_data = True)
ion.populate_ion_matrix()

ion.solve_no_source(np.array([1,0,0,0,0,0,0]),time)
ion.solve_time_independent()

plt.ion
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)

ax1.plot(time*1e3,ion.data['pops'][0,:,1,1],label='be0, met0',color='b')
ax1.hlines(ion.data['pops_ss'][0,0,1,1],0,10,color='b',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][1,:,1,1],label='be0, met1',color='g')
ax1.hlines(ion.data['pops_ss'][1,0,1,1],0,10,color='g',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][2,:,1,1],label='be1, met0',color='r')
ax1.hlines(ion.data['pops_ss'][2,0,1,1],0,10,color='r',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][3,:,1,1],label='be2, met0',color='c')
ax1.hlines(ion.data['pops_ss'][3,0,1,1],0,10,color='c',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][4,:,1,1],label='be2, met1',color='m')
ax1.hlines(ion.data['pops_ss'][4,0,1,1],0,10,color='m',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][5,:,1,1],label='be3, met0',color='y')
ax1.hlines(ion.data['pops_ss'][5,0,1,1],0,10,color='y',linestyle=':')

ax1.plot(time*1e3,ion.data['pops'][6,:,1,1],label='be4',color='k')
ax1.hlines(ion.data['pops_ss'][6,0,1,1],0,10,color='k',linestyle=':')

ax1.legend(fontsize='x-small')



#ion.solve_source(np.array([1,0,0,0,0,0,0]),np.array([1,0,0,0,0,0,0]),np.linspace(0,10,1e3))
