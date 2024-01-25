'''

Given cross-section data, this function calculates rate coefficients
given an electron energy distribution (EEDF)

cjperks, Jan 18th, 2024

TO DO:
    1) high-energy asymptotic cross-section for radiative
        recombination assumes you only have one ionized state

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
    DC_flag = False,
    limit = None,   # (optional), high-energy asymptotic behavior
    w_upr = None,   # (optional) Upper level total angular momentum
    w_lwr = None,   # (optional) Upper level total angular momentum
    ion_L = None,   # (optional) ionized state orbital angular momentum
    verbose = 1,
    use_rel = True, # flag to use relativistic corrections
    react = None,   # reaction type in ['ce', 'rr', 'ci']
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
            ntrans -- number of transitions considered
            NOTE: Assumed to be given on a sufficient grid so
                that setting any extrapolated values to zero
                causes negligible error

        engyXS -- [eV], dim(nE,ntrans), cross-section energy axis
            NOTE: Assumed to be not extremely fine, but fine enough
                so that a log-log interpolation causes negligible error

        m -- flag on definition of cross-section energy axis
            options: 
            1) if m=0, incident electron energy
            2) if m=1, scattered electron energy

        dE -- (optional), [eV], dim(ntrans,),
            if m=1 need the bound-bound transition energy difference

        DC_flag -- (optional), if XS is dielectronic capture strength
            Assumes cross-section is a delta function in energy
            NOTE: units of XS are then [eV *cm2]

    '''

    # Check
    if np.ndim(XS) == 1 and not DC_flag:
        XS = XS[:,None]
        engyXS = engyXS[:,None]

    # Prepares EEDF
    if EEDF == 'Maxwellian':
        EEDF, engyEEDF = _get_Max(
            Te=Te,
            use_rel=use_rel,
            ) # dim(ngrid, ntemp), ([1/eV], [eV])
    else:
        print('NON-MAXWELLIAN ELECTRONS NOT IMPLEMENTED YET')
        sys.exit(1)

    # If XS is dielectronic capture strength
    if DC_flag:
        return _calc_DC(
            EEDF=EEDF,
            engyEEDF=engyEEDF,
            XS=XS,
            engyXS=engyXS,
            m=m,
            dE=dE,
            use_rel=use_rel,
            ) # [cm3/s], dim(ntemp,ntrans)

    # Performs numerical integration
    else:
        return _calc_ratec(
            EEDF=EEDF,
            engyEEDF=engyEEDF,
            XS=XS,
            engyXS=engyXS,
            m=m,
            dE=dE,
            limit = limit,
            w_upr = w_upr,
            w_lwr = w_lwr,
            ion_L = ion_L,
            verbose = verbose,
            use_rel=use_rel,
            react=react,
            ) # [cm3/s], dim(ntemp,ntrans)


############################################################
#
#                       Utilities
#
############################################################

# Calculaes dielectronic capture rate coefficient
def _calc_DC(
    EEDF = None,        # [1/eV], dim(ngrid,ntemp)
    engyEEDF = None,    # [eV], dim(ngrid, ntemp)
    XS = None,          # [eV*cm2], dim(ntrans,)
    engyXS = None,      # [eV], dim(ntrans,)
    m = None,
    dE = None,          # [eV], dim(ntrans,)
    use_rel=None,
    ):

    # Init output
    ratec = np.zeros((EEDF.shape[1], len(XS))) # [cm3/s], dim(ntemp,ntrans)

    # Loop over transitions
    for nn in np.arange(len(XS)):
        if m == 0:
            engy_tmp = engyXS[nn]
        elif m == 1:
            engy_tmp = engyXS[nn] + dE[nn]

        # Incident electron velocity
        vel = _get_vel(
            E_inc = engyEEDF[:,tt],   
            use_rel = use_rel,
            ) # [cm/s], dim(ngrid,)

        # Loop over temperatures
        for tt in np.arange(EEDF.shape[1]):
            # Interpolates EEDF
            EEDF_tmp = 10**interp1d(
                np.log10(engyEEDF[:,tt]),
                np.log10(EEDF[:,tt]),
                bounds_error = False,
                fill_value = (-1e5,-1e5)
                )(np.log10(engy_tmp)) # [1/eV]

            # Calculates rate coefficient, [cm3/s]
            ratec[tt,nn] = XS[nn] *vel *EEDF_tmp

    # Output, [cm3/s], dim(ntemp,ntrans)
    return ratec

# Calculate rate coefficient integral
def _calc_ratec(
    EEDF = None,        # [1/eV], dim(ngrid,ntemp)
    engyEEDF = None,    # [eV], dim(ngrid, ntemp)
    XS = None,          # [cm2], dim(nE, ntrans)
    engyXS = None,      # [eV], dim(nE, ntrans), incident electron energy
    m = None,
    dE = None,          # [eV], dim(ntrans,)
    limit = None,       # dim(ntrans,2) if ce or dim(ntrans,4) if rr, ci
    w_upr = None,      
    w_lwr = None,
    ion_L = None,
    verbose = None,
    use_rel = None,
    react = None,
    ):

    # Init output
    ratec = np.zeros((EEDF.shape[1], XS.shape[1])) # [cm3/s], dim(ntemp,ntrans)
    if verbose == 1:
        XS_out = np.zeros(
            (EEDF.shape[0], EEDF.shape[1], XS.shape[1])
            ) # dim(ngrid, ntemp, ntrans)
        engy_out = np.zeros(
            (EEDF.shape[0], EEDF.shape[1])
            ) # dim(ngrid, ntemp)

    # Loop over temperatures
    for tt in np.arange(EEDF.shape[1]):
        # Incident electron velocity
        vel = _get_vel(
            E_inc = engyEEDF[:,tt],   
            use_rel = use_rel,
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
                fill_value = (-1e5,-1e5)
                )(np.log10(engyEEDF[:,tt])) # dim(ngrid,), [cm2]

            # Fill missing data right at threshold
            ## NOTE: collisional ioniz quickly drops to 0 near E_inc = dE
            ## NOTE: collisional excit is pretty nonlinear but flat enough this should be fine
            ## NOTE: radiative recomb blows up near E_inc = dE, so hopefully this is just a small error !!!
            if react != 'ci':
                indE = np.where(
                    (engyEEDF[:,tt] >= dE[nn])
                    & (engyEEDF[:,tt] <= engy_tmp[0])
                    )[0]
                XS_tmp[indE] = XS[0,nn]

            # Fill values with high-energy asymptotic behavior if available
            if limit is not None:
                XS_tmp = _get_limit(
                    XS_tmp = XS_tmp,
                    engy_tmp = engy_tmp,
                    engyEEDF = engyEEDF[:,tt],
                    dE = dE[nn],
                    limit = limit[nn,:],
                    w_upr = w_upr[nn],
                    w_lwr = w_lwr[nn],
                    ion_L = ion_L[0],
                    react = react,
                    )

            # Preforms integration
            ratec[tt,nn] = np.trapz(
                XS_tmp *vel * EEDF[:,tt],
                engyEEDF[:,tt]
                )

            if verbose == 1:
                XS_out[:,tt,nn] = XS_tmp
                engy_out[:,tt] = engyEEDF[:,tt] 

    # Output, [cm3/s], dim(ntemp,ntrans)
    if verbose == 1:
        return ratec, XS_out, engy_out
    else:
        return ratec

############################################################
#
#                       Calculation
#
############################################################

# Calculate high-energy asymptotic behavior of cross-section
## NOTE: So far, FAC-specific
def _get_limit(
    XS_tmp = None,
    engy_tmp = None,    # [eV], incident electron energy always
    engyEEDF = None,    # [eV], incident electron energy always
    dE = None,          # [eV], threshold energy
    limit = None,
    w_upr = None,
    w_lwr = None,
    ion_L = None,
    react = None,
    ):

    # Useful constants
    fine_struct2 = cnt.physical_constants['fine-structure constant'][0]**2 # []
    eV2Hartree = 1/cnt.physical_constants['Hartree energy in eV'][0] # [Hartree/eV]
    a02 = cnt.physical_constants['Bohr radius'][0]**2 *1e4 # [cm^2]

    # FAC cross-section data isn't really trustworthy about 10x transition energy
    tol = 10
    if engy_tmp[-1] >= tol*dE:
        indE = np.where(engyEEDF > tol*dE)[0]
    else:
        indE = np.where(engyEEDF> engy_tmp[-1])[0]

    # If unnecessary
    if len(indE) == 0:
        return XS_tmp
    
    if react == 'ce':
        # Electron kinetic momentum squared with fine structure correction
        k02 = 2 * engyEEDF[indE] *eV2Hartree *(
            1+
            0.5 *fine_struct2 *engyEEDF[indE] *eV2Hartree
            ) # [atomic units]

        # Bethe asymptotic form of collision strength for (first term) optically
        # allowed and (second term) forbidden transitions
        omega = (
            limit[0]
            * np.log(engyEEDF[indE]/dE)
            + limit[1]
            ) # []

        # Cross-section
        XS_tmp[indE] = (
            np.pi *omega
            / k02
            / (1 +2*w_lwr)
            ) * a02 # [cm2]

    elif react == 'rr':
        # Incident photon energy, [eV]
        E_gamma = (engyEEDF[indE] + dE)

        # Photon-electron energy, [eV]
        E_e = engyEEDF[indE]

        # Formula for bound-free oscillator strength, [1/Hartree]
        xx = (E_e + limit[3]) /limit[3]
        yy = (1 +limit[2]) /(np.sqrt(xx) + limit[2])

        dgf_dE = (
            E_gamma /(E_e +limit[3])
            * limit[0]
            * xx**(-3.5 -ion_L +0.5 *limit[1])
            * yy**limit[1]
            )

        eps = E_e *eV2Hartree # [Hartree]
        omega = E_gamma *eV2Hartree # [Hartree]

        # Function for photo-ionization cross-section, [atomic units]
        XS_PI = (
            2 *np.pi *np.sqrt(fine_struct2)
            /(1 +2*w_lwr)
            *(1 +fine_struct2 *eps)
            /(1 +0.5 *fine_struct2 *eps)
            *dgf_dE
            )

        # Function for radiative recombinarion cross-section, [cm2]
        XS_tmp[indE] = (
            fine_struct2/2
            * (1 +2*w_lwr)/(1 +2*w_upr)
            * omega**2
            / eps
            / (1 +0.5*fine_struct2 *eps)
            * XS_PI
            ) *a02

    elif react == 'ci':
        # Electron kinetic momentum squared with fine structure correction
        k02 = 2 * engyEEDF[indE] *eV2Hartree *(
            1+
            0.5 * fine_struct2 * engyEEDF[indE] *eV2Hartree
            ) # [atomic units]

        # Formula for collision strength
        xx = engyEEDF[indE]/dE
        yy = 1 - 1/xx

        omega = (
            limit[0] *np.log(xx)
            + limit[1] *yy**2
            + limit[2] *1/xx *yy
            + limit[3] *1/xx**2 *yy
            )

        # Cross-section
        XS_tmp[indE] = (
            omega
            / k02
            / (1 +2*w_lwr)
            ) * a02 # [cm2]

    # Output
    return XS_tmp

# Calculate Maxwellian energy distribution function
def _get_Max(
    Te=None, # [eV], dim(ntemp,)
    # Energy grid settings
    ngrid = int(1e3),    # energy grid resolution
    lwr = 1e-3,          # energy grid limits wrt multiple of Te
    upr = 5e1,           # !!!!!! Make sure these make sense for cross-sections
    use_rel = None,    # Use relativistic corrections
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

        # relativistic form
        if Te[tt] >= 10e3 and use_rel:
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

        # low-energy form
        else:
            EEDF[:,tt] = (
                2*np.sqrt(engy[:,tt]/np.pi)
                * Te[tt]**(-3/2)
                * np.exp(-engy[:,tt]/Te[tt])
                ) # dim(ngrid, ntemp); [1/eV] 

    # Output
    return EEDF, engy

# Incident electron velocity
def _get_vel(
    E_inc = None,   # [eV]
    use_rel = None
    ):

    # Useful constants
    eV2electron = cnt.e /(cnt.m_e*cnt.c**2) # [1/eV]

    # Relativistic form
    if use_rel:
        return (
            cnt.c *100 *np.sqrt(
                1 -
                1/(E_inc *eV2electron +1)**2
                ) 
            ) # [cm/s], dim(ngrid,)
    # Classical form
    else:         
        return (
            cnt.c *100 *np.sqrt(
                2*E_inc *eV2electron
                )
            ) # [cm/s], dim(ngrid,)
