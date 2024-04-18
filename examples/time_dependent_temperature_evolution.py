import os
import time
import numpy as np 
import matplotlib.pyplot as plt
from scipy import integrate
from colradpy.energy_balance import energy_balance


#%% Specify element to analyze and path to ACD and SCD adf11 files
element_name = 'carbon'
element_symbol = 'C'
num_charge_states = 7
year = '96'
pec_data_dir = r"C:\Users\Matt\Dropbox\Research\Auburn\Collisional radiative modelling\Atomic data"  # Need to specify path to ADAS ADF11 files on local machine
files = np.array([
    os.path.join(pec_data_dir, f'scd{year}_{element_symbol.lower()}.dat'),
    os.path.join(pec_data_dir, f'acd{year}_{element_symbol.lower()}.dat'),
])


#%% Specify input physics parameters for the energy balance
mass = 12.011 # mass of the species being modelled in amu
polarizability = 12.0 # dipole polarizabiliy in Bohr radii cubed of the modelled species' neutral state
electron_temp = 50 # background plasma electron temperature in eV
electron_dens = 1e13 # background plasma electron density in cm^-3
ion_temp = 75 # background plasma ion temperature in eV
ion_dens = 1e13 # background plasma ion density in cm^-3
ion_mass = 1.00794 # background plasma ion mass in amu
ion_charge_number = 1 # background plasma ion charge number
ion_species = "protium"


#%% Setup time-dependent energy balance
# Times for energy balance to be solved at. Unlike an ionization balance, the
# solution times for the energy balance need to begin sufficiently close to t=0
# that the populations and energies have not evolved significantly from their
# values at t=0. It helps to run an ionization balance first to see the time
# dynamics of the system and then use that info to refine the solution times
# for the energy balance.
t = np.geomspace(1e-10, 1e1, 301)
initial_abundances = np.zeros(num_charge_states)
initial_abundances[0] = 1 # 100% of the population in the neutral charge state
initial_abundances /= np.sum(initial_abundances)
initial_temperatures = np.zeros(num_charge_states)
initial_temperatures[0] = 1.5 # Neutral carbon sourced temperature characteristic of chemical erosion

balance = energy_balance(
    files,
    metas=np.array([0, 0, 0]),  # This argument isn't used when adf11 files are used
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

#%% Run time-dependent energy balance
t0 = time.perf_counter()
balance.solve()
t1 = time.perf_counter()
print(f"Time to calculate energy balance: {t1-t0:.1f} s")


#%% Calculate the time-averaged temperature for each species
average_energies = integrate.simpson(x=t, y=balance.data['processed']['energies'])
average_populations = integrate.simpson(x=t, y=balance.ion_balance.data['processed']['pops_td'][:, :, 0, 0])
average_temperatures = average_energies / (3/2 * average_populations)

# Remove bad data for charge states that aren't produced
bad_states = np.max(balance.ion_balance.data['processed']['pops_td'][:, :, 0, 0], axis=1) < 1e-6
average_energies[bad_states] = np.nan
average_temperatures[bad_states] = np.nan


#%% Plot time-dependent fractional abundances
fig, ax = plt.subplots(constrained_layout=True)
for charge_state in range(num_charge_states): # loop over all charge states
    ax.plot(
        t,
        balance.ion_balance.data['processed']['pops_td'][charge_state,:,0,0],
        label=f'$\mathrm{{{element_symbol}}}^{{{charge_state}\plus}}$',
    )
ax.set_xscale('log')
ax.set_ylim([0, 1])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Fractional abundance")
ax.set_title(
    f"Fractional ion balance for {element_name}\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  ",
)
leg = ax.legend()
leg.set_draggable(True)


#%% Plot time-dependent energies
fig, ax = plt.subplots(constrained_layout=True)
for charge_state in range(num_charge_states): # loop over all charge states
    ax.plot(
        balance.data['processed']['time'],
        balance.data['processed']['energies'][charge_state, :],
        label=f'$\mathrm{{{element_symbol}}}^{{{charge_state}\plus}}$',
    )
ax.set_xscale('log')
ax.set_ylim([0, None])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Energy (eV)")
ax.set_title(
    f"Energy balance for {element_name} in {ion_species} plasma\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  " + "\n"
    fr"$T_\mathrm{{i}}$ = {balance.data['user']['ion_temp']:.1f} eV "
    fr"$n_\mathrm{{i}}$ = {balance.data['user']['ion_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
)
leg = ax.legend()
leg.set_draggable(True)


#%% Plot time-dependent temperatures
population_threshold = 1e-5  # 1e-6 is limit where you start seeing numerical instabilities

fig, ax = plt.subplots(constrained_layout=True)
for charge_state in range(num_charge_states): # loop over all charge states
    temperature = balance.data['processed']['temperatures'][charge_state, :].copy()
    is_valid = balance.ion_balance.data['processed']['pops_td'][charge_state,:,0,0] < population_threshold
    temperature[is_valid] = np.nan
    ax.plot(
        balance.data['processed']['time'],
        temperature,
        label=f'$\mathrm{{{element_symbol}}}^{{{charge_state}\plus}}$',
    )
ax.set_xscale('log')
ax.set_ylim([0, 1.05*ion_temp])
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (eV)")
ax.set_title(
    f"Temperature evolution for {element_name} in {ion_species} plasma\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  " + "\n"
    fr"$T_\mathrm{{i}}$ = {balance.data['user']['ion_temp']:.1f} eV "
    fr"$n_\mathrm{{i}}$ = {balance.data['user']['ion_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
)
leg = ax.legend()
leg.set_draggable(True)

#%% Plot emissivity-weighted average temperatures
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(np.arange(num_charge_states), average_temperatures, marker="o")
ax.set_xlim([-1/2, num_charge_states-1/2])
ax.set_xlabel("Charge state")
ax.set_ylabel("Temperature (eV)")
ax.set_title(
    f"Emissivity-averaged temperatures for {element_name} in {ion_species} plasma\n"
    fr"$T_\mathrm{{e}}$ = {balance.data['user']['electron_temp']:.1f} eV  "
    fr"$n_\mathrm{{e}}$ = {balance.data['user']['electron_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  " + "\n"
    fr"$T_\mathrm{{i}}$ = {balance.data['user']['ion_temp']:.1f} eV  "
    fr"$n_\mathrm{{i}}$ = {balance.data['user']['ion_dens']:5.1e} $\mathrm{{cm}}^{{-3}}$  "
)
