import numpy as np
from colradpy.r8necip import *
####################################################################################################
def ecip_rates(energy,ion_pot,zpla,zpla1,charge_state,temp_grid):
    """This function gives energy levels to r8necip for ECIP ionization to be
calculated for and then returns ECIP rates for each energy level on the temperature
grid provided.

    :param energy: Energy levels
    :type energy: array

    :param ion_pot: ionization potential
    :type ion_pot: float

    :param zpla: ionization potential
    :type zpla: array

    :param zpla1: ionization potential
    :type zpla1: array

    :param charge_state: charge state
    :type charge_state: int

    :param temp_grid: Temperature grid to return ECIP rates on
    :type temp_grid: array

    :returns: array --ECIP rates on the temperature grid

    """
    ecip = np.zeros((len(energy),len(ion_pot),len(temp_grid)))
    for p in range(0,len(ion_pot)):
        #there was a stupid problem beacuse {X} only designates ecip not other ionization        
        ion_inds = np.where( zpla[:,p] > -1)[0]
        ion_inds2 = np.where( zpla1 == p +1)[0]
        ion_inds = np.intersect1d(ion_inds,ion_inds2)

        ecip[ion_inds,p,:] = r8ecip(charge_state, ion_pot[p],
                         energy[ion_inds],zpla[ion_inds,p],
                          temp_grid*11604.5)
    return ecip
