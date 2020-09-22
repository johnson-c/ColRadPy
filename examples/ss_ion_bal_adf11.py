import numpy as np 
import matplotlib.pyplot as plt
path_to_colradpy = '/home/curtis/git/tmp_tmp/'
import sys
sys.path.append(path_to_colradpy)
from colradpy import ionization_balance
#loading in the prerequisit packages need numpy as matplotlib always
#need to append the path of where colradpy is located on local machine
#we will be running an ionization balance so this needs to be imported



#need to specify the location of the SCD and ACD adf11 files on the local machine
fils = np.array(['/home/curtis/imladris/gcr_for_class/scd89_fe.dat',
                 '/home/curtis/imladris/gcr_for_class/acd89_fe.dat'])


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







