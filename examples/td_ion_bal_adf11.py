import numpy as np 
import matplotlib.pyplot as plt
from os.path import exists
from os import mkdir
from pathlib import Path
import sys
from urllib import request
# Import ColRadPy
sys.path.append(
    str(Path(__file__).parent.parent)
)
from colradpy import ionization_balance


# Variables
# Set up output folders
EXAMPLES_PATH: Path = Path(__file__).parent
EXAMPLES_INPUT_PATH: Path = EXAMPLES_PATH / "input"
EXAMPLES_OUTPUT_PATH: Path = EXAMPLES_PATH / "output"
OUTPUT_PATH: Path = EXAMPLES_OUTPUT_PATH / Path(__file__).name.split('.')[0]
ADAS_PATH: Path = EXAMPLES_INPUT_PATH / "Open-ADAS"
INPUT_PATH: Path = ADAS_PATH / "adas#26"
# Making output directories
paths = [EXAMPLES_OUTPUT_PATH, OUTPUT_PATH, ADAS_PATH, INPUT_PATH]
for p in paths:
    if not exists(p):
        mkdir(p)
#loading in the prerequisit packages need numpy as matplotlib always
#need to append the path of where colradpy is located on local machine
#we will be running an ionization balance so this needs to be imported



#need to specify the location of the SCD and ACD adf11 files on the local machine
adf11_files = [
    'scd89_fe.dat',
    'acd89_fe.dat'
]
adf11_urls = [
    'https://open.adas.ac.uk/download/adf11/scd89/scd89_fe.dat',
    'https://open.adas.ac.uk/download/adf11/acd89/acd89_fe.dat'
]

for adf11_file, url in zip(adf11_files, adf11_urls):
    adf11_path: Path = INPUT_PATH / adf11_file
    if not exists(adf11_path):
        with request.urlopen(url) as f:
            adf11: str = f.read().decode('utf-8')
            with open(adf11_path, "x") as f_adf11:
                f_adf11.write(adf11)

fils = np.array([INPUT_PATH / f for f in adf11_files], dtype=str)

#####################################
# Time dependent ionization balance #
#####################################

temp_grid = np.array([8])#temperature for the TD ion balance
dens_grid = np.array([6.31e11])# density for the TD ion balance

time = np.geomspace(1.e-6,2,1000000) #times for the ionization to be done at
                                    #note that because the matrix exponentiation method is used
                                    #any time step or number of times can be chosen
initial_abundances = np.zeros(27) # intitial abundances for the ionization balance
initial_abundances[0] = 1.        #Here 100% of the population states in the neutral charge state

ion_td = ionization_balance(fils, 
                            np.array([0,0,0]),
                            temp_grid,dens_grid,
                            adf11_files=True,
                            soln_times=time, #solution times to be calculated
                            init_abund = initial_abundances) #the initial abundance at t=0

ion_td.populate_ion_matrix() #
ion_td.solve_no_source()




###############################
# Plotting for time dependent #
###############################

fig2 = plt.figure(figsize=(12,6), facecolor = "white")
for i in range(0,5):
    plt.plot(time,ion_td.data['processed']['pops_td'][i,:,0,0])

plt.text(ion_td.data['user']['soln_times'][-1] * 0.08, 0.5, r"Final Time = " +
         "{:5.3f}".format(ion_td.data['user']['soln_times'][-1]) +' S', fontsize = 16)
plt.text(ion_td.data['user']['soln_times'][-1] * 0.08, 0.4, r"$N_e$ = " + 
         "{:5.3e}".format(ion_td.data['user']['dens_grid'][0]) + r"  $cm^{-3}$", fontsize = 16)

plt.text(ion_td.data['user']['soln_times'][-1] * 0.08, 0.3, r"$T_e$ = 8 eV" , fontsize = 16)


plt.xscale('log')
plt.xlabel("Time (s)", fontsize = 16)
plt.ylabel("Fractional Abundance", fontsize = 16)
plt.title("Time Dependent Fractional Ion Balance for Iron", fontsize = 16)

fig2.savefig(OUTPUT_PATH / f"{OUTPUT_PATH.name}.pdf", format='pdf')
