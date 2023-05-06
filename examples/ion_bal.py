#######################################################################
#
# Note that this script requires the C adf04 datafiles from openADAS
#
#######################################################################
#
# Imports
import numpy as np
import matplotlib.pyplot as plt
from os.path import exists
from os import mkdir
from pathlib import Path
import sys
from urllib import request
# Import ColRadPy
sys.path.append('/home/curtis/git/tmp_tmp/')
from colradpy import ionization_balance


# Variables
# Set up output folders
EXAMPLES_PATH: Path = Path(__file__).parent
EXAMPLES_INPUT_PATH: Path = EXAMPLES_PATH / "input"
EXAMPLES_OUTPUT_PATH: Path = EXAMPLES_PATH / "output"
OUTPUT_PATH: Path = EXAMPLES_OUTPUT_PATH / Path(__file__).name.split('.')[0]
ADAS_PATH: Path = EXAMPLES_INPUT_PATH / "Open-ADAS"
INPUT_PATH: Path = ADAS_PATH / "adas#6"
# Making output directories
paths = [EXAMPLES_OUTPUT_PATH, OUTPUT_PATH, ADAS_PATH, INPUT_PATH]
for p in paths:
    if not exists(p):
        mkdir(p)
# Plotting
plt.ion()
plt.rc('font',size=8)
plt.rcParams['font.weight'] = 'semibold'
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


# Download ADF04 from Open-ADAS
adf04_files = []
charges = range(6)
for charge in charges:
    url: str = f"https://open.adas.ac.uk/download/adf04/adas][6/mom97_ls][c{charge}.dat"
    adf04_file: str = f"mom97_ls#c{charge}.dat"
    adf04_files.append(adf04_file)
    adf04_path: Path = INPUT_PATH / adf04_file
    if not exists(adf04_path):
        with request.urlopen(url) as f:
            adf04: str = f.read().decode('utf-8')
            with open(adf04_path, "x") as f_adf04:
                f_adf04.write(adf04)

fils = np.array([INPUT_PATH / f for f in adf04_files], dtype=str)


temp = np.array([100])
dens= np.array([1.e13])
htemp = np.array([10])
hdens = np.array([1.e13])

time = np.geomspace(1.e-6,200,100000)

metas = [np.array([0]),np.array([0,1]),np.array([0,1]),
         np.array([0]),np.array([0,1]),np.array([0]),np.array([0])]#number of metastable
                           #this should match the metastables at the top of the adf04 file
                           #this information is used to calculate the QCD values
                           #without it only the SCD, ACD and XCD for a species will be calculated

initial_abundances = np.array([1,0,0,0,0,0,0,0,0,0]) #initial fractional abundances of the charge states
                                               #everything is starting the in neutral ground state
                                               #this array must be the same size as the total
                                               #number of metastables in the system
source = np.array([0,0,0,0,0,0,0,0,0,0])

ion_0 = ionization_balance(fils, metas, temp,
                         dens,
                         keep_charge_state_data = False,
                         use_cx=False,
                               temp_dens_pair=False,scale_file_ioniz=True,
                           soln_times = time, init_abund = initial_abundances,
                           source=source)

ion_0.populate_ion_matrix()
ion_0.solve_time_independent()
ion_0.solve_no_source()


source = np.array([1,0,0,0,0,0,0,0,0,0])
ion_1 = ionization_balance(fils, metas, temp,
                         dens,
                         keep_charge_state_data = False,
                         use_cx=False,
                               temp_dens_pair=False,scale_file_ioniz=True,
                           soln_times = time, init_abund = initial_abundances,
                           source=source)

ion_1.populate_ion_matrix()
ion_1.solve_time_independent()
ion_1.solve_source()




########################################################################################################################
initial_abundances = np.array([1e13,0,0,0,0,0,0,0,0,0]) #initial fractional abundances of the charge states
                                               #everything is starting the in neutral ground state
                                               #this array must be the same size as the total
                                               #number of metastables in the system
source = np.array([1e13,0,0,0,0,0,0,0,0,0])
ion_e = ionization_balance(fils, metas, temp,
                         dens,
                         keep_charge_state_data = False,
                         use_cx=False,
                               temp_dens_pair=False,scale_file_ioniz=True,
                           soln_times = time, init_abund = initial_abundances,
                           source=source)

ion_e.populate_ion_matrix()
ion_e.solve_time_independent()
ion_e.solve_source()



initial_abundances = np.array([0,0,0,0,0,0,0,0,0,0]) #initial fractional abundances of the charge states
                                               #everything is starting the in neutral ground state
                                               #this array must be the same size as the total
                                               #number of metastables in the system
source = np.array([1,0,0,0,0,0,0,0,0,0])
ion_s = ionization_balance(fils, metas, temp,
                         dens,
                         keep_charge_state_data = False,
                         use_cx=False,
                               temp_dens_pair=False,scale_file_ioniz=True,
                           soln_times = time, init_abund = initial_abundances,
                           source=source)

ion_s.populate_ion_matrix()
ion_s.solve_time_independent()
ion_s.solve_source()


#This plot shows that the electron density does not impact the time
#required to reach equilbrium. The time required to reach equilbrium
#is only dependent on the ratio of the initial abundances to the source term

fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.93,left=0.11,right=0.99)
ca = np.array(['b','r','g','y','c','olivedrab','deeppink','purple','lime','chocolate'])

for i in range(0,len(initial_abundances)):

    plt.plot(time, ion_1.data['processed']['pops_td_source'][i,:,0,0]/\
             np.sum(ion_1.data['processed']['pops_td_source'][:,:,0,0],axis=0),color=ca[i],linestyle='--',label=str(i))
    plt.plot(time, ion_e.data['processed']['pops_td_source'][i,:,0,0]/\
             np.sum(ion_e.data['processed']['pops_td_source'][:,:,0,0],axis=0),color='k',linestyle=':')
    
ax1.set_ylabel('Ionization Fraction (-)',weight='semibold')
ax1.set_xlabel('Time (s)',weight='semibold')
plt.title('T$_e$ = 100 eV  n$_e$ = 10$^{13}$ cm$^{-1}$',weight='semibold')
plt.legend()
plt.xscale('log')

fig.savefig(OUTPUT_PATH / "carbon_ion_bal.pdf", format='pdf')
    


#This plot shows the difference between the two extreme cases
#solid line is no source term. The dashed line has only a source term.
#source terms are into the ground state only for this example
fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
fig.subplots_adjust(bottom=0.15,top=0.93,left=0.11,right=0.99)

for i in range(0,len(initial_abundances)):
    plt.hlines( ion_0.data['processed']['pops_ss'][i,0,0]/\
                np.sum(ion_0.data['processed']['pops_ss'][:,0,0],axis=0),0,time[-1],color=ca[i])
    plt.plot(time, ion_s.data['processed']['pops_td_source'][i,:,0,0]/\
             np.sum(ion_s.data['processed']['pops_td_source'][:,:,0,0],axis=0),color=ca[i],linestyle='--')

ax1.set_ylabel('Ionization Fraction (-)',weight='semibold')
ax1.set_xlabel('Time (s)',weight='semibold')
plt.title('T$_e$ = 100 eV  n$_e$ = 10$^{13}$ cm$^{-1}$',weight='semibold')
plt.legend()
plt.xscale('log')
