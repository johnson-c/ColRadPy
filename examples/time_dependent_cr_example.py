import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
from os import mkdir
from pathlib import Path
import sys
# Import ColRadPy
sys.path.append(
    str(Path(__file__).parent.parent)
)
from colradpy import colradpy


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

#Time dependent CR modeling
td_t = np.geomspace(1.e-5,.01,1000)
td_n0 = np.zeros(30)
td_n0[0] = 1.

fil = str(EXAMPLES_INPUT_PATH / 'cpb03_ls#be0.dat') #adf04 file
temperature_arr = np.array([10]) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e13]) # cm-3
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
title = 'Time dependent solution of CR Be I no source term'
ax1.set_title(title)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Population (-)')

file_name = title.replace(' ', '-')
fig.savefig(OUTPUT_PATH / f"{file_name}.pdf", format='pdf')





td_t = np.geomspace(1.e-5,2,1000)
td_n0 = np.zeros(30)
td_n0[0] = 1.
td_s = np.zeros(30)
td_s[0] = 1.



be_s = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
              use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True,
              td_t=td_t,td_n0=td_n0,td_source=td_s)

be_s.solve_cr()
be_s.solve_time_dependent()

fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.92,left=0.115,right=0.965)
plt.plot(be_s.data['user']['td_t'],
         be_s.data['processed']['td']['td_pop'][0,:,0,0],
         label='Ground')
plt.plot(be_s.data['user']['td_t'],
         be_s.data['processed']['td']['td_pop'][1,:,0,0],
         label='level 1')
plt.plot(be_s.data['user']['td_t'],
         be_s.data['processed']['td']['td_pop'][-1,:,0,0],
         label='ion')
ax1.legend(fontsize='x-small',loc='best')
title = 'Time dependent solution of CR Be I with source term'
ax1.set_title(title)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Population (-)')


file_name = title.replace(' ', '-')
fig.savefig(OUTPUT_PATH / f"{file_name}.pdf", format='pdf')


