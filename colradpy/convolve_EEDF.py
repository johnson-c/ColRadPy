'''

Given cross-section data, this function calculates rate coefficients
given an electron energy distribution (EEDF)

cjperks, Jan 18th, 2024


'''

# Modules
import numpy as np
import scipy.constants as cnt
from scipy.special import kn
from scipy.interpolate import interp1d

############################################################
#
#                       Main
#
############################################################

def convolve_EEDF(
    EEDF = None,    # [1/eV]
    Te = None,      # [eV], dim(ntemp,)
    XS = None,      # [cm2], dim(nE,ntrans)
    engyXS = None,  # [eV], dim (nE,ntrans)
    m = 0,          # 0 or 1
    dE = 0,         # [eV]
    ):
    '''
    INPUTS:
        EEDF -- [1/eV], 
            options:
            1) 'Maxwellian' --> pre-defined maximal entropy distribution function
                if Te <10 keV = Maxwell-Boltzmann
                if Te >10 keV = Maxwell-Juttner (relativistic)

            2) NOT IMPLEMENTED YET (non-Maxwellians)

        Te -- [eV], dim(ntemp,), characteristic temperature of EEDF
            if using non-Maxwellian then whatever characteristic 
            dimension that isn't the energy axis

        XS -- [cm2], dim(nE,ntrans), cross-section data
            nE -- number of energy grid points
            ntrnas -- number of transitions considered

        engyXS -- [eV], dim(nE,ntrans), cross-section energy axis
            NOTE: Assumed to be not extremely fine

        m -- flag on definition of cross-section energy axis
            options: 
            1) if m=0, incident electron energy
            2) if m=1, scattered electron energy

        dE -- (optional), [eV], dim(ntrans,),
            if m=1 need the bound-bound transition energy difference

    '''

    # Useful constants
    eV2electron = cnt.e /(cnt.m_e*cnt.c**2) # [1/eV]

    # Check
    if ndim(XS) == 1:
        XS = XS[:,None]
        engyXS = engyXS[:,None]

    # Prepares EEDF
    if EEDF == 'Maxwellian':
        EEDF, engyEEDF = _get_Max(
            Te=Te,
            ) # dim(ngrid, ntemp), ([1/eV], [eV])
    else:
        print('NON-MAXWELLIAN ELECTRONS NOT IMPLEMENTED YET')
        sys.exit(1)

    # Init output
    ratec = np.zeros((len(Te))) # [cm3/s], dim(ntemp,ntrans)

    # Loop over temperatures
    for tt in np.arange(len(Te)):
        # Incident electron velocity, relativistic form
        vel =  cnt.c *100 *np.sqrt(
            1 -
            1/(engyEEDF[:,tt] *eV2electron +1)**2
            ) # [cm/s], dim(ngrid,)

        # Loop over transitions
        for nn in np.arange(XS.shape[1]):
            # Interpolates cross-section onto finer grid of EEDF
            if m == 0:
                engy_tmp = engyXS[:,nn]
            elif m == 1:
                engy_tmp = engyXS[:,nn] + dE[nn]

            XS_tmp = 10**interp1d(
                np.log10(engy_tmp),
                np.log10(XS[:,nn]),
                bounds_error=False,
                fill_value = (0,0)
                )(np.log10(engyEEDF[:,tt])) # dim(ngrid,), [cm2]

            # Preforms integration
            ratec[tt,nn] = np.trapz(
                XS_tmp *vel * EEDF[:,tt],
                engyEEDF[:,tt]
                )

    # Output, [cm3/s], dim(ntemp,ntrans)
    return ratec


############################################################
#
#                       Utilities
#
############################################################

# Calculate Maxwellian energy distribution function
def _get_Max(
    Te=None, # [eV], dim(ntemp,)
    # Energy grid settings
    ngrid = int(1e4),    # energy grid resolution
    lwr = 1e-3,          # energy grid limits wrt multiple of Te
    upr = 5e1,           # !!!!!! Make sure these make sense for cross-sections
    ):
    '''
    NOTE: EEDF energy axis defined as incident electron energy
    '''

    # Useful constants
    eV2electron = cnt.e /(cnt.m_e*cnt.c**2) # [1/eV]

    # Init output
    EEDF = np.zeros((ngrid, len(Te))) # dim(ngrid, ntemp), [1/eV]
    engy = np.zeros((ngrid, len(Te))) # dim(ngrid, ntemp), [eV]

    # Loop over temperature
    for tt in np.arange(len(Te)):
        engy[:,tt] = np.logspace(
            np.log10(Te[tt]*lwr), 
            np.log10(Te[tt]*upr),
            ngrid
            ) # [eV]

        # low-energy form
        if Te[tt] <= 10e3:
            EEDF[:,tt] = (
                2*np.sqrt(Einc/np.pi)
                * temps[tt]**(-3/2)
                * np.exp(-Einc/temps[tt])
                ) # dim(ngrid, ntemp); [1/eV] 

        # relativistic form
        else:
            # Factors
            theta = Te[tt] * eV2electron
            gamma = 1 + engy[:,tt] * eV2electron
            beta = np.sqrt(1-1/gamma**2)

            EEDF[:,tt] = (
                gamma**2
                * beta
                /theta
                /kn(2, 1/theta)
                *np.exp(-gamma/theta)
                ) * eV2electron # dim(ngrid, ntemp); [1/eV]

    # Output
    return EEDF, engy

