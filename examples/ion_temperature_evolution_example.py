import numpy as np 
import matplotlib.pyplot as plt
from colradpy.energy_balance import energy_balance


#%% Specify element to analyze and path to ACD and SCD adf11 files
path_to_atomic_data = ""  # Path on local machine where adf11 files are stored
element_name = 'carbon'
element_symbol = 'C'
num_charge_states = 7
year = '96'
files = np.array([
    f'{path_to_atomic_data}/scd{year}_{element_symbol.lower()}.dat',
    f'{path_to_atomic_data}/acd{year}_{element_symbol.lower()}.dat',
])


#%% Specify input physics parameters for the energy balance
dens = 1.2e11 # total density in cm^-3 of the species being modelled summed over all charge states
mass = 12.011 # mass of the species being modelled in amu
polarizability = 12.0 # dipole polarizabiliy in Bohr radii cubed of the modelled species' neutral state
electron_temp = 100 # background plasma electron temperature in eV
electron_dens = 1e13 # background plasma electron density in cm^-3
ion_temp = 100 # background plasma ion temperature in eV
ion_dens = 1e13 # background plasma ion density in cm^-3
ion_mass = 1.007 # background plasma ion mass in amu
ion_charge_number = 1 # background plasma ion charge number
ion_species = "protium"


#%% Run time-dependent ionization balance
# Times for energy balance to be recorded at. Note that the solver picks its own time steps internally.
t = np.geomspace(1e-8, 1e-0, 301)
initial_abundances = np.zeros(num_charge_states)
initial_abundances[0] = 1 # 100% of the population in the neutral charge state
initial_temperatures = np.zeros(num_charge_states)
initial_temperatures[0] = 1.5 # Neutral carbon sourced at temperature characteristic of chemical erosion

balance = energy_balance(
    files,
    metas=np.array([0, 0, 0]),  # This argument isn't used when adf11 files are used
    dens=dens,
    mass=mass,
    polarizability=polarizability,
    electron_temp=electron_temp,
    electron_dens=electron_dens,
    ion_temp=ion_temp,
    ion_dens=ion_dens,
    ion_mass=ion_mass,
    ion_charge_number=ion_charge_number,
    adf11_files=True,
    soln_times=t,
    init_abund=initial_abundances,
    init_temp=initial_temperatures,
)
balance.solve()


#%% Plot time-dependent fractional abundances
fig, ax = plt.subplots(constrained_layout=True)
for i in range(num_charge_states): # loop over all charge states
    ax.plot(
        t,
        balance.ion_balance.data['processed']['pops_td'][i,:,0,0],
        label=f'{element_symbol} {i}+',
    )
ax.set_xscale('log')
ax.set_ylim([0, 1])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Fractional abundance")
ax.set_title(
    f"Fractional ion balance for {element_name}\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
)
leg = ax.legend()
leg.set_draggable(True)


#%% Plot time-dependent energies
fig, ax = plt.subplots(constrained_layout=True)
for i in range(num_charge_states): # loop over all charge states
    ax.plot(
        balance.data['processed']['time'],
        balance.data['processed']['energies'][i, :],
        label=f'{element_symbol} {i}+',
    )
ax.set_xscale('log')
ax.set_ylim([0, None])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Energy (eV)")
ax.set_title(
    f"Energy balance for {element_name}\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$T_\mathrm{{i}}$ = {balance.data['user']['ion_temp']:.1f} eV" + "\n"
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
    fr"$n_\mathrm{{i}}$ = {balance.data['user']['ion_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
    fr"$n_\mathrm{{{element_symbol}}}$ = {balance.data['user']['dens']:5.1e} $\mathrm{{cm}}^{{-3}}$"
)
leg = ax.legend()
leg.set_draggable(True)


#%% Plot time-dependent temperatures
fig, ax = plt.subplots(constrained_layout=True)
for i in range(num_charge_states): # loop over all charge states
    ax.plot(
        balance.data['processed']['time'],
        balance.data['processed']['temperatures'][i, :],
        label=f'{element_symbol} {i}+',
    )
ax.set_xscale('log')
ax.set_ylim([0, 1.05*ion_temp])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (eV)")
ax.set_title(
    f"Temperature evolution for {element_name}\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$T_\mathrm{{i}}$ = {balance.data['user']['ion_temp']:.1f} eV" + "\n"
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
    fr"$n_\mathrm{{i}}$ = {balance.data['user']['ion_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
    fr"$n_\mathrm{{{element_symbol}}}$ = {balance.data['user']['dens']:5.1e} $\mathrm{{cm}}^{{-3}}$"
)
leg = ax.legend()
leg.set_draggable(True)
