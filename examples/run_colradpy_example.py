import sys
sys.path.append('../') #starting in 'examples' so need to go up one
from colradpy import colradpy
import numpy as np
from os.path import exists
from os import mkdir
from pathlib import Path


# Variables
# Set up output folders
EXAMPLES_PATH: Path = Path(__file__).parent
EXAMPLES_INPUT_PATH: Path = EXAMPLES_PATH / "input"
EXAMPLES_OUTPUT_PATH: Path = EXAMPLES_PATH / "output"
OUTPUT_PATH: Path = EXAMPLES_OUTPUT_PATH / Path(__file__).name.split('.')[0]
# Making output directories
paths = [EXAMPLES_OUTPUT_PATH, OUTPUT_PATH]
for p in paths:
    if not exists(p):
        mkdir(p)

fil = str(EXAMPLES_INPUT_PATH / 'cpb03_ls#be0.dat') #adf04 file
temperature_arr = np.linspace(1,100,100) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e13,4.e14]) # cm-3

#calling the colradpy class with the various inputs
be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True, 
              use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)

be.solve_cr() #solve the CR equations with the quasistatic method


#the block of code below is what 'solve_cr()' is doing
if(be.data['user']['use_ionization']):
    be.make_ioniz_from_reduced_ionizrates()
if(be.data['user']['suppliment_with_ecip']):
    be.make_ecip()
    be.suppliment_with_ecip()
if(be.data['user']['use_recombination']):
    be.make_recombination_rates_from_file()
if(be.data['user']['use_recombination_three_body']):
    be.make_three_body_recombination()
be.make_electron_excitation_rates()
be.populate_cr_matrix()
be.solve_quasi_static()
#end of solve_cr()












'''



import matplotlib.pyplot as plt
plt.ion()
met = 0 #metastable 0, this corresponds to the ground state
te = 0 #first temperature in the grid
ne = 0 #frist density in the grid

fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.105,right=0.965)
ax1.vlines(be.data['processed']['wave_air'],
           np.zeros_like(be.data['processed']['wave_air']),
           be.data['processed']['pecs'][:,met,te,ne])
ax1.set_xlim(0,1000)
ax1.set_title('PEC spectrum  T$_e$=' +str(be.data['user']['temp_grid'][te])+\
              ' eV  ne=' + "%0.*e"%(2,be.data['user']['dens_grid'][ne]) + ' cm$^{-3}$',size=10)
ax1.set_xlabel('Wavelength (nm)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')


print(np.shape(be.data['processed']['wave_air']),
      np.shape(be.data['processed']['pecs']),
      np.shape(be.data['processed']['pec_levels']))
#(320,) (320, 3, 1, 1) (320, 2)

upper_ind = 7 #ninth excited state
lower_ind = 0  #ground state

pec_ind = np.where( (be.data['processed']['pec_levels'][:,0] == upper_ind) &\
                    (be.data['processed']['pec_levels'][:,1] == lower_ind))[0]


#plot the temeprature dependence of the chosen pec at first density in the grid
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.93,left=0.105,right=0.965)
ax1.set_title('Temperature dependence of line ' +\
              str(be.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
ax1.plot(be.data['user']['temp_grid'],be.data['processed']['pecs'][pec_ind[0],met,:,ne])
ax1.set_xlabel('Temperature (eV)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')

#plot the density dependence of the chosen pec at first density in the grid
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.93,left=0.125,right=0.965)
ax1.set_title('Density dependence of line ' +\
              str(be.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
ax1.plot(be.data['user']['dens_grid'],be.data['processed']['pecs'][pec_ind[0],met,te,:])
ax1.set_xlabel('Density (cm$^{-3}$)')
ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')


#want to find the index of Be I line at 351.55
pec_ind = np.where( (be.data['processed']['wave_air'] <352) &\
                    (be.data['processed']['wave_air'] >351))
print('Wavelength from file ' + str(be.data['processed']['wave_air'][pec_ind[0]]))
#Wavelength from file [351.55028742]
print('PEC upper and lower levels '+ str(be.data['processed']['pec_levels'][pec_ind[0]]))
#PEC upper and lower levels [[25  2]]
'''
