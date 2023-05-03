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


######################################################
# Time independent ionization balance (steady state) #
######################################################


temp_grid = np.geomspace(1,40000,1000) #temperature grid in (eV) for steady state ionization balance
dens_grid = np.array([6.31e11]) # density grid in (cm-3) for steady state ionization balance

#call ionization balance class the second required input needs to be left as a blank array
#this is for other functionality that when running from adf04 files insteady of adf11
ion = ionization_balance(fils, #a text array of the files used for the ionization balance
                         np.array([]),#leave blank when using adf11 files like this example
                         temp_grid, #temperature grid to do the ionization balance over
                         dens_grid, #density grid to do the ionization balance over
                         adf11_files=True) #flag needs to be true when using adf11 files as input

ion.populate_ion_matrix() #take the rates interpolated from the adf11 file and put then in the ionization balance matrix
ion.solve_time_independent()#solve the ionization balance matrix in steady state



#############################
# Plotting for steady state #
#############################

fig1 = plt.figure(figsize=(12, 6), facecolor='white')
for i in range(0,np.shape(ion.data['ion_matrix'])[0]):#plot all of the charge states in the ionization balance
    plt.semilogx(ion.data['user']['temp_grid'], ion.data['processed']['pops_ss'][i,:,0])
    plt.xlabel("Temperature (eV)", fontsize = 16)
    plt.ylabel("Fractional Abundance", fontsize = 16)
    plt.title("Equilibrium Fractional Ion Balance for Iron", fontsize = 16)


    plt.text(1e3, 0.9, r"$N_e$ = " + "{:5.3e}".format(ion.data['user']['dens_grid'][0]),
             fontsize = 16)#print the density of the ionization to the plot

plt.axvline(np.array([8]), color="black", linestyle="--")#a time dependent ionization balance will be run at this temperature
plt.text(8 * 1.1, 0.9, r"$T_e$ for time dependent analysis")

fig1.savefig(OUTPUT_PATH / f"{OUTPUT_PATH.name}.pdf", format='pdf')
