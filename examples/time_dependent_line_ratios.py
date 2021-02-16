import sys
sys.path.append('../')
from colradpy import colradpy
import numpy as np
import matplotlib.pyplot as plt


#############################################
# Setting time dependent CR modeling params #
#############################################
td_t = np.geomspace(1.e-10,1.e-3,100000)
td_n0 = np.zeros(30) 
td_n0[0] = 1.#only poulation starts in the ground state

#####################################
# Setting non-time-dependent params #
#####################################
fil = 'cpb03_ls#be0.dat' #adf04 file
temperature_arr = np.array([10]) #eV
metastable_levels = np.array([0])   #metastable level, just ground chosen here
density_arr =     np.array([1.e13]) # cm-3





###################
# Call to ColRadPy#
###################

#no recombination here to make things easier for normalization but can be done too
be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=False,
              use_recombination_three_body = False,use_ionization=True,suppliment_with_ecip=True,
              td_t=td_t,td_n0=td_n0)
be.solve_cr() #solve the QS
be.solve_time_dependent() #solve the TD




###################################
# Comparing TD and QS populations #
###################################
qs_ratio_pop = be.data['processed']['pops'][0,:,0,0]/be.data['processed']['pops'][1,:,0,0]

#must integrate over the total time eplased here just do it for a long time so it will be
# the same as equilbrium. If non-steady state choose some smaller time
td_ratio_pop = np.trapz(be.data['processed']['td']['td_pop'][1,:,0,0],be.data['user']['td_t'])/\
         np.trapz(be.data['processed']['td']['td_pop'][2,:,0,0],be.data['user']['td_t'])

print('Quasistatic population level_1/level_2:  ',qs_ratio_pop )
print('Time dependent population integration over time level_1/level_2:  ',td_ratio_pop)


############################
# Comparing TD and QS PECs #
############################

#The ratio of PECs can also be taken, just take the ratio of the two largest PECs
# PEC 1 and PEC 33
#standard ratio of popuations
qs_ratio_pecs = be.data['processed']['pecs'][1,:,0,0]/be.data['processed']['pecs'][33,:,0,0]

#Must integrate over the total time eplased here just do it for a long time so it will be
#the same as equilbrium. If non-steady state choose some smaller time
td_ratio_pecs = np.trapz(be.data['processed']['td']['pecs'][1,:,0,0],be.data['user']['td_t'])/\
                np.trapz(be.data['processed']['td']['pecs'][33,:,0,0],be.data['user']['td_t'])  


print('Quasistatic PEC ratio PEC_1/PEC_33:  ',qs_ratio_pecs )
print('Time dependent PEC integration over time PEC_1/PEC_33:  ',td_ratio_pecs)


#####################################
# Plotting for time dependent levels#
#####################################

plt.ion()

plt.rc('font',size=8)
plt.rcParams['font.weight'] = 'semibold'
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)



fig, ax1 = plt.subplots(1,1,figsize = (6.5,4.5),dpi=300)
fig.subplots_adjust(bottom=0.11,top=0.94,left=0.18,right=0.97,hspace=0.3,wspace=.35)
ax1.tick_params(axis='both',direction='in')

ax1.plot(be.data['user']['td_t'],be.data['processed']['td']['td_pop'][0,:,0,0],label='Ground')
ax1.plot(be.data['user']['td_t'],be.data['processed']['td']['td_pop'][1,:,0,0],label='Metastable')
ax1.plot(be.data['user']['td_t'],be.data['processed']['td']['td_pop'][-1,:,0,0],label='Ionized')
plt.xscale('log')
ax1.set_xlabel('Time (s)',weight='semibold')
ax1.set_ylabel('Fractional population (-)',weight='semibold')
plt.legend()
plt.tight_layout()
