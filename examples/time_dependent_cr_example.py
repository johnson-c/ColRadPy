import sys
sys.path.append('../')
from colradpy_class import colradpy
import numpy as np
import matplotlib.pyplot as plt

#Time dependent CR modeling
td_t = np.geomspace(1.e-5,.1,1000)
td_n0 = np.zeros(30)
td_n0[0] = 1.

fil = 'cpb03_ls#be0.dat' #adf04 file
temperature_arr = np.array([10]) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e9]) # cm-3
be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
              use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True,
              td_t=td_t,td_n0=td_n0)
be.solve_cr()
be.solve_time_dependent()

fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.1,right=0.965)
plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][0,:,0,0],
         label='Ground')
plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][1,:,0,0],
         label='level 1')

plt.plot(be.data['user']['td_t'],
         (be.data['processed']['td']['td_pop'][2,:,0,0]/ be.data['processed']['td']['td_pop'][0,:,0,0])/np.sum((be.data['processed']['td']['td_pop'][2,:,0,0]/ be.data['processed']['td']['td_pop'][0,:,0,0])),
         label='level 1')




plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][-1,:,0,0],
         label='ion')
ax1.legend(fontsize='x-small',loc='best')
ax1.set_title('Time dependent solution of CR Be I no source term')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Population (-)')






td_t = np.geomspace(1.e-5,1,1000)
td_n0 = np.zeros(30)
td_n0[0] = 1.
td_s = np.zeros(30)
td_s[0] = 1.
fil = 'cpb03_ls#be0.dat' #adf04 file
temperature_arr = np.array([10]) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e8]) # cm-3
be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
              use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True,
              td_t=td_t,td_n0=td_n0,td_source=td_s)

be.solve_cr()
be.solve_time_dependent()

fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.115,right=0.965)
plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][0,:,0,0],
         label='Ground')
plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][1,:,0,0],
         label='level 1')
plt.plot(be.data['user']['td_t'],
         be.data['processed']['td']['td_pop'][-1,:,0,0],
         label='ion')
ax1.legend(fontsize='x-small',loc='best')
ax1.set_title('Time dependent solution of CR Be I with source term')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Population (-)')
