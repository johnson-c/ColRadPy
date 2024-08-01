from functools import partial
import numpy as np
from scipy import constants
from scipy.integrate import solve_ivp
from colradpy.ionization_balance_class import ionization_balance
from colradpy.solve_matrix_exponential import eval_matrix_exponential_solution


class energy_balance(object):
    """
    Solves the time-dependent evolution of the temperatures of a species charge
    states as it is ionized in a plasma. This problem is essentially the same
    as an ionization balance, with the addition of heating of the charge states
    by collisions with the main ion species in the plasma.

    This class makes use of 'ionization_balance_class.py' to first perform an
    ionization balance calculation before performing the energy balance
    calculation. Unlike an ionization balance calculation, the energy balance
    involves solving a nonlinear system of rate equations, due to the nonlinear
    temperature dependence of the collisional energy transfer between plasma
    species. So, the efficient analytic matrix ODE solver in
    'solve_matrix_exponential.py' cannot be used. The energy balance can
    therefore only be performed for a single set of main ion/electron
    densities/temperatures at a time (it is not vectorized across a grid of
    temperatures/densities like the ionization balance).

    :Args:       
      :param fils: Array of the input adf04 or adf11 files
      :type fils: string array

      :param metas: List of arrays for the metastable levels in a charge state
      :type metas: list

      :param mass: Mass of the species being modelled in (amu)
      :type mass: float
      
      :param polarizability: Dipole polarizability of the neutral species being modelled in Bohr radii cubed
      :type polarizability: float

      :param electron_temp: Electron temperature of the background plasma in (eV)
      :type electron_temp: float

      :param electron_dens: Electron density of the background plasma in (cm^-3)
      :type electron_dens: float

      :param ion_temp: Main ion temperature of the background plasma in (eV)
      :type ion_temp: float

      :param ion_dens: Main ion density of the background plasma in (cm^-3)
      :type ion_dens: float

      :param ion_mass: Main ion species mass in (amu)
      :type ion_mass: float

      :param ion_charge_number: Main ion species charge number
      :type ion_charge_number: int
      
      :param init_temp: Initial temperatures in (eV) of the charge states for the modelled species
      :type init_temp: float array

      :param kwargs: Keyword arguments to pass to `ionization_balance`
    """

    def __init__(self, fils, metas, mass, polarizability, electron_temp,
                 electron_dens, ion_temp, ion_dens, ion_mass, ion_charge_number,
                 init_temp=np.array([]), **kwargs):

        # Set up ionization balance 
        self.ion_balance = ionization_balance(
            fils, metas, temp_grid=np.array([electron_temp]),
            dens_grid=np.array([electron_dens]), **kwargs,
        )

        # Initialize data
        self.data = {}
        self.data["user"] = {}
        self.data["user"]["mass"] = mass
        self.data["user"]["polarizability"] = polarizability
        self.data["user"]["electron_temp"] = electron_temp
        self.data["user"]["electron_dens"] = electron_dens
        self.data["user"]["ion_temp"] = ion_temp
        self.data["user"]["ion_dens"] = ion_dens
        self.data["user"]["ion_mass"] = ion_mass
        self.data["user"]["ion_charge_number"] = ion_charge_number
        self.data["user"]["init_temp"] = init_temp


    def solve(self, n0=np.array([]), T0=np.array([]), td_t=np.array([])):
        """
        Solves the time-dependent energy balance given initial populations and temperatures

        :param n0: The initial fractional abundance of each charge state at t=0
        :type n0: float array
        
        :param T0: The intial temperature in (eV) of each charge state at t=0
        :type T0: float array

        :param td_t: Solution times
        :type td_t: float array

        populates the abundances, energies, and temperatures for all solution times
        """
        # Set default parameter values
        if(n0.size < 1):
            n0 = self.ion_balance.data['user']['init_abund']
        if(T0.size < 1):
            T0 = self.data['user']['init_temp']
        if(td_t.size < 1):
            td_t = self.ion_balance.data['user']['soln_times']

        # Solve the ionization balance to get fractional abundances vs time
        self.ion_balance.populate_ion_matrix()
        self.ion_balance.solve_no_source(n0, td_t)

        # Calculate initial total energy of each charge state
        E0 = 3/2 * n0 * T0

        # Set up the energy matrix to track energy transfer due to ionization/recombination processes
        self.energy_matrix = (
            self.ion_balance.data["ion_matrix"][:, :, 0, 0]
            * self.data["user"]["electron_dens"]
        )

        # Calculate Debye length (assume mobility of electrons is higher than all other species)
        self.debye_length = debye_length(
            charges=np.array([-1]),
            densities=np.array([self.data["user"]["electron_dens"]]),
            temperatures=np.array([self.data["user"]["electron_temp"]]),
        )

        # Create function that returns population in each charge state at any time
        pops_fun = partial(
            eval_matrix_exponential_solution,
            n0=n0,
            eigenvalues=self.ion_balance.data["processed"]["eigen_val"],
            eigenvectors=self.ion_balance.data["processed"]["eigen_vec"],
        )

        # Calculate solution using scipy solver
        sol = solve_ivp(
            self._evolve_energy, t_span=[0, td_t[-1]], t_eval=td_t, y0=E0,
            args=(pops_fun,), method="BDF",
        )

        # Package solution data
        self.data["processed"] = {}
        self.data["processed"]["time"] = td_t
        self.data["processed"]["abundances"] = self.ion_balance.data["processed"]["pops_td"][:, :, 0, 0]
        self.data["processed"]["energies"] = sol.y
        self.data["processed"]["temperatures"] = self.data["processed"]["energies"] / (3/2 * self.data["processed"]["abundances"])
        self.data["processed"]["nfev"] = sol.nfev
        self.data["processed"]["status"] = sol.status
        self.data["processed"]["message"] = sol.message
        self.data["processed"]["success"] = sol.success


    def _evolve_energy(self, time, energies, pops_fun):
        """
        Calculates the time derivative of the energy in each charge state.

        :param time: Time at which to compute the energy time derivative
        :type time: float
        
        :param energies: Energy in each charge state at the current time
        :type energies: float array
        
        :param pops_fun: Function that returns the population in each charge state for an input time
        :type pops_fun: fun(t)

        Returns a vector of the time derivative of the energy in each charge state
        """
        # Get fractional abundance for this time and calculate state temperatures
        pops = pops_fun(np.array([time]))[:, 0, 0, 0]
        # Handle special cases where populations are zero (typically at beginning)
        temps = np.divide(energies, 3/2 * pops, out=np.zeros_like(energies), where=(pops != 0))
        temps[temps < 0] = 0  # Prevent numerical noise from making negative temperatures

        # Set up energy transfer matrix & source vector (starting with ionization matrix)
        energy_matrix = self.energy_matrix.copy()
        num_charge_states = energy_matrix.shape[0]
        energy_source = np.zeros(num_charge_states)

        # Create vector of energy transfer frequencies for each charge state
        frequencies = np.zeros(num_charge_states)
        frequencies[0] = neutral_energy_transfer_frequency(
            np.array([self.data["user"]["mass"], self.data["user"]["ion_mass"]]) * constants.value("atomic mass constant"),
            self.data["user"]["ion_dens"] * 1e6,  # cm^-3 -> m^-3
            self.data["user"]["polarizability"],
        )
        frequencies[1:] = plasma_energy_transfer_frequency(
            masses=np.array(
                [self.data["user"]["mass"], self.data["user"]["ion_mass"]]
            ) * constants.value("atomic mass constant"),  # amu -> kg
            charges=np.column_stack((
                np.arange(1, num_charge_states),
                np.full(num_charge_states - 1, self.data["user"]["ion_charge_number"]),
            )),
            densities=np.full(num_charge_states - 1, self.data["user"]["ion_dens"] * 1e6),  # cm^-3 -> m^-3
            temperatures=np.column_stack((
                temps[1:],
                np.full(num_charge_states - 1, self.data["user"]["ion_temp"]),
            )),
            debye_length=self.debye_length,
        )
        frequencies[~np.isfinite(frequencies)] = 0 # Handle special cases where temperature is zero

        # Collisional energy exchange is proportional to the temperature difference
        # between the two species: dE/dt = -frequency * density (temperature_diff).
        # To allow reuse of the ionization matrix, the negative temperature term
        # (representing energy loss from the species due to collisions with main
        # ions) is included in the energy matrix, while the positive temperature
        # term (representing heating of the species due to collisions with main
        # ions) is included in an energy source vector that is then added to the
        # product of the energy matrix with the energies vector.
        energy_matrix[np.diag_indices(num_charge_states)] -= 2/3 * frequencies
        energy_source = frequencies * pops * self.data["user"]["ion_temp"]
        energy_derivative = np.einsum("ij,j->i", energy_matrix, energies) + energy_source

        return energy_derivative


def plasma_energy_transfer_frequency(
    masses, charges, densities, temperatures, debye_length
):
    """
    Calculates the collisional energy transfer frequency between two plasma
    species.
    
    The thermal equilibration of two species is given by:
        3/2 * n_1 * dT_1/dt = -nu_1,2 * n_1 * (T_1 - T_2)
    where n is density, T is temperature, t is time, and nu is the collisional
    energy transfer frequency. This frequency characterizes the rate at which
    energy is transferred by collisions from species 2 to species 1. It is
    given by Eq. 2.4.55 in Jim Callen's Plasma Kinetic Theory lecture notes.
    Note that the frequency given in the more widely known NRL Plasma Formulary
    differs from the definition here by a factor of 3/2 (because that frequency
    gives the temperature relaxation rate, while this formula gives the energy
    transfer rate). To convert this frequency into the temperature relaxation
    frequency, multiply the result by 2/3.

    Parameters
    ----------
    masses : array, shape (..., 2)
        Mass of each species in kilograms. 
    charges : array, shape (..., 2)
        Charge number of each species.
    densities : array, shape (...)
        Density of species 2 in particles per cubic meter.
    temperatures : array, shape (..., 2)
        Temperature of each species in electron-volts.
    debye_length : array, shape (...)
        Debye length of the plasma in meters.

    Returns
    -------
    frequency : array
        Frequency in inverse seconds characterizing the timescale for which
        energy is collisionally transferred from species 2 to species 1.
    """
    thermal_velocities = np.sqrt(2 * constants.e * temperatures / masses)
    log = coulomb_logarithm(masses, charges, temperatures, debye_length)
    frequency = (
        (
            densities * charges[..., 0]**2 * charges[..., 1]**2
            * constants.e**4 * log
        )
        / (
            constants.pi**(3/2) * constants.epsilon_0**2 * masses[..., 0]
            * masses[..., 1] * ((thermal_velocities**2).sum(axis=-1))**(3/2)
        )
    )
    return frequency


def neutral_energy_transfer_frequency(masses, plasma_density, neutral_dipole_polarizability):
    """
    Calculates the collisional energy transfer frequency between a neutral and
    a plasma species.

    This frequency characterizes the rate at which energy is transferred by
    collisions from a plasma species to a neutral species. It is given by Eq. 9
    in [J.D. Hey et al., J. Phys. B: At. Mol. Opt. Phys. 35, 1525-1553 (2002)].

    Note that this formula is an elastic scattering approximation that does not
    include contributions from longer range interactions.

    The thermal equilibration of the neutral species with the plasma species is
    given by:
        3/2 * n_n * dT_n/dt = -nu_n,p * n_n * (T_n - T_p)
    where n_n is neutral density, T_n/T_p is neutral/plasma temperature, t is
    time, and nu is the collisional energy transfer frequency.

    Parameters
    ----------
    masses : array, shape (..., 2)
        Mass of the neutral and plasma species in kilograms. 
    plasma_density : array, shape (...)
        Density of the plasma species in particles per cubic meter.
    neutral_dipole_polarizability : array, shape (...)
        Dipole polarizability of the neutral species in Bohr radii cubed.
        Polarizabilities vary between metastable states and can be found for
        many atomis in [V.P. Shevelko, Atoms and their Spectroscopic
        Properties, Berlin: Springer (1997)].

    Returns
    -------
    frequency : array
        Frequency in inverse seconds characterizing the timescale for which
        energy is collisionally transferred from the plasma to neutral species.
    """
    reduced_mass = masses[..., 0] * masses[..., 1] / (masses[..., 0] + masses[..., 1])
    momentum_transfer_rate_coefficient = (
        2.21 * constants.pi * constants.alpha * constants.c
        * np.sqrt(
            constants.m_e / reduced_mass
            * neutral_dipole_polarizability * constants.value("Bohr radius")**4
        )
    )
    frequency = (
        2 * plasma_density * reduced_mass / (masses[..., 0] + masses[..., 1])
        * momentum_transfer_rate_coefficient
    )
    return 3/2 * frequency  # Include 3/2 factor to go from temperature to energy


def coulomb_logarithm(masses, charges, temperatures, debye_length):
    """
    Calculates the coulomb logarithm for collisions between two plasma species.

    This equation comes from Eqs. 2.1.19 and 2.2.21-23 in Jim Callen's Plasma
    Kinetic Theory lecture notes.
    
    Parameters
    ----------
    masses : array, shape (..., 2)
        Mass of each species in kilograms.
    charges : array, shape (..., 2)
        Charge number of each species.
    temperatures : array, shape (..., 2)
        Temperature of each species in electron-volts.
    debye_length : array, shape (...)
        Debye length of the plasma in meters.

    Returns
    -------
    coulomb_logarithm : array
    """
    reduced_mass = masses[..., 0] * masses[..., 1] / masses.sum(axis=-1)
    thermal_velocities = np.sqrt(2 * constants.e * temperatures / masses)
    min_impact_parameter_classical = (
        np.abs(charges[..., 0] * charges[..., 1]) * constants.e**2
        / (
            4 * np.pi * constants.epsilon_0
            * reduced_mass * (thermal_velocities**2).sum(axis=-1) / 2
        )
    )
    min_impact_parameter_quantum = (
        constants.hbar
        / (reduced_mass * np.sqrt((thermal_velocities**2).sum(axis=-1)))
    )
    min_impact_parameter = np.max(
        [min_impact_parameter_classical, min_impact_parameter_quantum],
        axis=0,
    )
    return np.log(debye_length / min_impact_parameter)


def debye_length(charges, densities, temperatures):
    """
    Calculates the Debye length of a multi-species plasma.

    charges : array, shape (..., S)
        Charge number of each species. Species axis must be last.
    densities : array, shape (..., S)
        Density of each species in particles per cubic meter. Species axis
        must be last.
    temperatures : array, shape (..., S)
        Temperature of each species in electron-volts. Species axis must be
        last.

    Returns
    -------
    debye_length : array, shape (...)
        Debye length of the plasma in meters.
    """
    #TODO: The Debye length calculation
    #should not include plasma species for which the temperature is so low that
    #the q * phi / k_B * T << 1 approximation becomes invalid. Not entirely
    #sure how to check this approximation... I think it would involve an
    #iterative process where the analytic solution for the shielded potential
    #with all species included is used to evaluate q_s * phi / T_s at
    #r = (n_s)^(-1/3). If for any species this expression is not << 1, then
    #that species should be removed, the Debye length recalculated, and so on
    #until all species involved in the calculation satisfy this requirement.
    return (
        (densities * charges**2 * constants.e**2)
        / (constants.epsilon_0 * constants.e * temperatures)
    ).sum(axis=-1)**(-1/2)


# Tests for collisional energy transfer rates
if __name__ == "__main__":
    proton_to_carbon2_rate = plasma_energy_transfer_frequency(  # Energy transfer rate from protons to C2+
        masses=np.array([12.011 * constants.value("atomic mass constant"), constants.proton_mass]),
        charges=np.array([2, 1]),
        densities=np.array([1e17, 5e18]),
        temperatures=np.array([1.5, 15]),
    )
    electron_to_carbon2_rate = plasma_energy_transfer_frequency(  # Energy transfer rate from electrons to C2+
        masses=np.array([12.011 * constants.value("atomic mass constant"), constants.electron_mass]),
        charges=np.array([2, -1]),
        densities=np.array([1e17, 5e18]),
        temperatures=np.array([1.5, 15]),
    )
    proton_to_carbon0_rate = neutral_energy_transfer_frequency(
        masses=np.array([12.011 * constants.value("atomic mass constant"), constants.proton_mass]),
        plasma_density=5e18,
        neutral_dipole_polarizability=12.0,
    )
    electron_to_carbon0_rate = neutral_energy_transfer_frequency(
        masses=np.array([12.011 * constants.value("atomic mass constant"), constants.electron_mass]),
        plasma_density=5e18,
        neutral_dipole_polarizability=12.0,
    )
    deuteron_to_electron_rate = plasma_energy_transfer_frequency(
        masses=np.array([constants.electron_mass, constants.value("deuteron mass")]),
        charges=np.array([-1, 1]),
        densities=np.array([1e21, 1e21]),
        temperatures=np.array([1e3, 10e3]),
    )
