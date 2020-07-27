import sys
sys.path.append('../')
import numpy as np
import matplotlib.pyplot as plt
from colradpy import ionization_balance

fils = np.array(['cpb03_ls#be0.dat','cpb03_ls#be1.dat','be2_adf04','be3_adf04']) #the adf04 files
temp = np.linspace(5,100,5) #temp grid
dens = np.array([1.e11,1.e14]) #density grid
metas = [np.array([0,1]),np.array([0]),np.array([0,1]),np.array([0])]#number of metastable
                           #this should match the metastables at the top of the adf04 file
                           #this information is used to calculate the QCD values
                           #without it only the SCD, ACD and XCD for a species will be calculated

time = np.linspace(0,.01,1.e4) #solution time grid

initial_abundances = np.array([1,0,0,0,0,0,0]) #initial fractional abundances of the charge states
                                               #everything is starting the in neutral ground state
                                               #this array must be the same size as the total
                                               #number of metastables in the system




ion = ionization_balance(fils, metas, temp, dens,
                         soln_times = time, init_abund = initial_abundances,
                         keep_charge_state_data = True,use_cx=False) #ionization balance call

ion.populate_ion_matrix() #populating the matrix with the GCR coefficients
                                               
ion.solve_no_source() #solving the system with no source of particles

ion.solve_time_independent() #solve the time independent equations to get the limits
                             #for the time dependent case


plt.ion
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.90,left=0.11,right=0.99)

ax1.plot(time*1e3,ion.data['processed']['pops_td'][0,:,1,1],label='be0, met0',color='b')#time dependent
ax1.hlines(ion.data['processed']['pops_ss'][0,1,1],0,10,color='b',linestyle=':')#time indpenent limit

ax1.plot(time*1e3,ion.data['processed']['pops_td'][1,:,1,1],label='be0, met1',color='g')
ax1.hlines(ion.data['processed']['pops_ss'][1,1,1],0,10,color='g',linestyle=':')

ax1.plot(time*1e3,ion.data['processed']['pops_td'][2,:,1,1],label='be1, met0',color='r')
ax1.hlines(ion.data['processed']['pops_ss'][2,1,1],0,10,color='r',linestyle=':')

ax1.plot(time*1e3,ion.data['processed']['pops_td'][3,:,1,1],label='be2, met0',color='c')
ax1.hlines(ion.data['processed']['pops_ss'][3,1,1],0,10,color='c',linestyle=':')

ax1.plot(time*1e3,ion.data['processed']['pops_td'][4,:,1,1],label='be2, met1',color='m')
ax1.hlines(ion.data['processed']['pops_ss'][4,1,1],0,10,color='m',linestyle=':')

ax1.plot(time*1e3,ion.data['processed']['pops_td'][5,:,1,1],label='be3, met0',color='y')
ax1.hlines(ion.data['processed']['pops_ss'][5,1,1],0,10,color='y',linestyle=':')

ax1.plot(time*1e3,ion.data['processed']['pops_td'][6,:,1,1],label='be4',color='k')
ax1.hlines(ion.data['processed']['pops_ss'][6,1,1],0,10,color='k',linestyle=':')

ax1.legend(fontsize='x-small')

ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Fractional Abundance (-)')
plt.title('Temperature 28 eV, Density 1*10$^{14}$ cm$^{-3}$')










temp = np.linspace(2,100,200) #temp grid
ion = ionization_balance(fils, metas, temp, dens, keep_charge_state_data = False,use_cx=False)
ion.populate_ion_matrix()
ion.solve_time_independent()


fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.99,left=0.11,right=0.99)

ax1.plot(temp,ion.data['processed']['pops_ss'][0,:,1],label='be0, met0',color='b')

ax1.plot(temp,ion.data['processed']['pops_ss'][1,:,1],label='be0, met1',color='g')

ax1.plot(temp,ion.data['processed']['pops_ss'][2,:,1],label='be1, met0',color='r')

ax1.plot(temp,ion.data['processed']['pops_ss'][3,:,1],label='be2, met0',color='c')

ax1.plot(temp,ion.data['processed']['pops_ss'][4,:,1],label='be2, met1',color='m')

ax1.plot(temp,ion.data['processed']['pops_ss'][5,:,1],label='be3, met0',color='y')

ax1.plot(temp,ion.data['processed']['pops_ss'][6,:,1],label='be4',color='k')

ax1.legend(fontsize='x-small')

ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('Fractional Abundance (-)')






ion.solve_source(np.array([1,0,0,0,0,0,0]),np.array([1,0,0,0,0,0,0]),np.linspace(0,.06,1e3))
