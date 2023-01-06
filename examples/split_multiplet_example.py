import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
from os import mkdir
from pathlib import Path
import sys
from urllib import request
# Import ColRadPy
sys.path.append('../')
from colradpy import colradpy


# Variables
# Set up output folders
EXAMPLES_PATH: Path = Path(__file__).parent
EXAMPLES_INPUT_PATH: Path = EXAMPLES_PATH / "input"
EXAMPLES_OUTPUT_PATH: Path = EXAMPLES_PATH / "output"
OUTPUT_PATH: Path = EXAMPLES_OUTPUT_PATH / Path(__file__).name.split('.')[0]
ADAS_PATH: Path = EXAMPLES_INPUT_PATH / "Open-ADAS"
INPUT_PATH: Path = ADAS_PATH / "adas#2"
# Making output directories
paths = [EXAMPLES_OUTPUT_PATH, OUTPUT_PATH, ADAS_PATH, INPUT_PATH]
for p in paths:
    if not exists(p):
        mkdir(p)


# Download ADF04 from Open-ADAS
url: str = f"https://open.adas.ac.uk/download/adf04/adas][2/mom97_ls][he0.dat"
adf04_file: str = 'mom97_ls#he0.dat'
adf04_path: Path = INPUT_PATH / adf04_file
if not exists(adf04_path):
    with request.urlopen(url) as f:
        adf04: str = f.read().decode('utf-8')
        with open(adf04_path, "x") as f_adf04:
            f_adf04.write(adf04)

fil = str(adf04_path)  #adf04 file location

temp_grid = np.array([500]) #temperature grid in (eV) for steady state ionization balance
dens_grid = np.array([1.e13]) # density grid in (cm-3) for steady state ionization balance
meta = np.array([0]) #index of the chosen metastables here only the ground is chosen

he = colradpy(fil,  #fil location
               meta, #metastable indexes
               temp_grid, #temperature grid
               dens_grid, #density grid
               use_ionization=True, #ionization rates included
               suppliment_with_ecip=True, #ECIP filling in where ionization not available
               use_recombination_three_body=False, #Do not use 3body recombination
               use_recombination=False,   # Do not use recombination
               rate_interp_col='cubic'
              )

he.make_electron_excitation_rates() #interpolate eciation rates on the user defined grid
he.make_ioniz_from_reduced_ionizrates() #interpolate ionization rates on the user defined grid
he.suppliment_with_ecip() #make ECIP rates where ionization not available
he.populate_cr_matrix()  #place rates into the CR matrix
he.solve_quasi_static()  #solve the CR matrix

he.get_nist_levels_txt()
he.split_structure_terms_to_levels()
he.shift_j_res_energy_to_nist()
he.split_pec_multiplet()



#Set up pretty plot settings
plt.rc('font',size=8)
plt.rcParams['font.weight'] = 'semibold'
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


fig, ax1 = plt.subplots(1,1,figsize = (6.5,3),dpi=300)
fig.subplots_adjust(bottom=0.08,top=0.94,left=0.08,right=0.97,hspace=0.3,wspace=.35)
ax1.tick_params(axis='both',direction='in')




ax1.vlines(he.data['processed']['wave_air'],
           np.zeros_like(he.data['processed']['wave_air']),
           he.data['processed']['pecs'][:,0,0,0],label='Term resolved')

ax1.vlines(he.data['processed']['split']['wave_air'],
           np.zeros_like(he.data['processed']['split']['wave_air']),
           he.data['processed']['split']['pecs'][:,0,0,0],label='Level resolved',color='r')

ax1.set_xlabel('Wavelength (nm)',weight='semibold')
ax1.set_ylabel('Intensity (AU)',weight='semibold')
ax1.set_xlim(587.5,587.7)
plt.ylim(0,1e-10)
plt.legend()
plt.tight_layout()
fig.savefig(OUTPUT_PATH / "split_multiplet_example_mom97_ls#he0.pdf", format='pdf')
