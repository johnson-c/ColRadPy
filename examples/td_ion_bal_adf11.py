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
