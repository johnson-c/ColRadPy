import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/curtis/git/tmp_tmp/')
from colradpy import ionization_balance
#######################################################################
#
#Note that this script requires the C adf04 datafiles from openADAS
#
#######################################################################

plt.ion()
plt.rc('font',size=8)
plt.rcParams['font.weight'] = 'semibold'
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


path = './adas#6/'
fils = np.array([path+'mom97_ls#c0.dat',path+'mom97_ls#c1.dat',path+'mom97_ls#c2.dat',
                 path+'mom97_ls#c3.dat',path+'mom97_ls#c4.dat',path+'mom97_ls#c5.dat'],dtype=str)


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
