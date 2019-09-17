import sys
sys.path.append('../')
from colradpy import colradpy
import numpy as np

fil = 'cpb03_ls#be1.dat' #adf04 file
temperature_arr = np.linspace(1,100,20) #eV
metastable_levels = np.array([0,1])   #ground and level 1 chosen to be metastable
density_arr =     np.array([1.e13,8.e13,4.e14]) # cm-3
beii = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
              use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)
beii.solve_cr()

#plotting the QCD
import matplotlib.pyplot as plt
plt.ion
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.125,right=0.965)
ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['qcd'][0,1,:,0]*1e5,
         label = 'metastable cross coupling coefficient 1->2')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['qcd'][1,0,:,0]*1e5,
         label = 'metastable cross coupling coefficient 2->1')
ax1.legend()
ax1.set_title('QCD plot')
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('QCD * 10$^5$ (cm$^{-3}$ s$^{-1}$)')

#plotting the SCD
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.125,right=0.965)
ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['scd'][0,0,:,0],
         label = 'metastable cross coupling coefficient 1->1+')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['scd'][0,1,:,0],
         label = 'metastable cross coupling coefficient 1->2+')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['scd'][1,0,:,0],
         label = 'metastable cross coupling coefficient 2->1+')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['scd'][1,1,:,0],
         label = 'metastable cross coupling coefficient 2->2+')

ax1.legend(fontsize='x-small',loc='best')
ax1.set_title('SCD plot')
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('SCD (ion cm$^{-3}$ s$^{-1}$)')

#plotting the ACD
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.075,right=0.965)
ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['acd'][0,0,:,0],
         label = 'metastable cross coupling coefficient 1+->1')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['acd'][0,1,:,0],
         label = 'metastable cross coupling coefficient 2+->1')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['acd'][1,0,:,0],
         label = 'metastable cross coupling coefficient 1+->2')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['acd'][1,1,:,0],
         label = 'metastable cross coupling coefficient 2+->2')

ax1.legend(fontsize='x-small',loc='best')
ax1.set_title('ACD plot')
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('ACD (rec cm$^{-3}$ s$^{-1}$)')

#plotting the XCD
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.12,right=0.965)
ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['xcd'][0,1,:,0],
         label = 'metastable cross coupling coefficient 1+->2+')

ax1.plot(beii.data['user']['temp_grid'],
         beii.data['processed']['scd'][1,0,:,0],
         label = 'metastable cross coupling coefficient 2+->1+')
ax1.legend(fontsize='x-small',loc='best')
ax1.set_title('XCD plot')
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('XCD (cm$^{-3}$ s$^{-1}$)')
