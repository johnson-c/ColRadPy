'''

read_FAC.py is a subfunction that reads atomic data created by
the Flexible Atomic Code (FAC) to be used in collisional-
radiative modeling in ColRadPy

NOTE: Assumes you have a local version of FAC

cjperks
Dec 18, 2023

TO DO:
    1) Add non-Maxwellian convolution
    2) Missing energy level data (zpla, zpla1)
    3) Dielectronic recombination/autoionization data
    4) Charge exchange data

'''

# Module
from pfac import rfac, crm
import os
import numpy as np
import copy
import colradpy.read_FAC_utils as utils

############################################################
#
#                       Main
#
############################################################

def read_FAC(
    # Ion species of interest because not stored in FAC data files....
    ele = None,         # Species name
    nele = None,        # Number of electrons
    Zele = None,        # Nuclear charge
    # File Management
    fil = None,         # Common path to FAC files, excluding physics extentions
    # Physics controls
    EEDF = None,        # if None -> assumes Maxwell-averages rates from pfac.fac.MaxwellRate
    reacts = None,      # if None -> looks for all file suffixes
    Te = None,          # if not using MaxwellRate files, [eV], dim(ntemp,)
    verbose = 1,
    ):

    ######## -------- Determines which files to search for -------- ########

    # Electron energy distribution settings
    if EEDF == 'Maxwellian_mr':
        # Use Maxwell-averaged rates from pfac.fac.MaxwellRate
        use_mr = True
    elif EEDF == 'Maxwellian':
        use_mr = False
    else:
        print('NON-MAXWELLIAN ELECTRONS NOT IMPLEMENTED YET')
        sys.exit(1)

    # FAC data files to search for
    if reacts == 'incl_all':
        reacts = [
            'en',      # Energy levels
            'tr',      # Einstein coefficients
            'ce',      # Collisional excitation
            'rr',      # Radiative recombination
            #'ai',      # Autoionization/dielectronic recombination
            'ci',      # Collision ionization
            ]

    ######## -------- Reads data -------- ########

    # Initialize output
    FAC = {}
    FAC['rates'] = {}

    # Energy levels
    if 'en' in reacts:
        FAC = _en(
            FAC=FAC,
            fil=fil,
            ele=ele,
            nele=nele,
            Zele=Zele,
            )
    # Error check
    else:
        print('NEED TO INCLUDE ENERGY LEVEL DATA IN MODELING!!!')
        sys.exit(1)

    # Gets rate coefficient data
    # Includes: collisional excit, radiative recomb, collisional ioniz
    if use_mr:
        # Use Maxwell-averaged rates from pfac.fac.MaxwellRate
        FAC = _get_mr(
            FAC=FAC,
            fil=fil,
            verbose=verbose,
            reacts=reacts,
            )
    else:
        # Use cross-section data files
        FAC = _get_xs(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            Te=Te,
            verbose=verbose,
            reacts=reacts,
            )

    # Einstein coefficients
    if 'tr' in reacts:
        FAC = _tr(
            FAC=FAC,
            fil=fil,
            nele=nele,
            Zele=Zele,
            )
    # Error check
    else:
        print('NEED TO INCLUDE EINSTEIN COEFFICIENT DATA IN MODELING!!!')
        sys.exit(1)

    # Autoionization/dielectronic recombination
    if 'ai' in reacts:
        FAC = _ai(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            Te=Te,
            )

    # Charge exchange
    if 'cx' in reacts:
        print('CHARGE EXCHANGE IS NOT IMPLEMENTED YET!!!')
    # If empty
    else:
        FAC['rates']['cx'] = {}
        FAC['rates']['cx']['cx_transitions'] = np.asarray([])
        FAC['rates']['cx']['cx_excit'] = np.asarray([])

    # Infinite energy points
    # I assume you made your FAC data setting the appropriate energies
    FAC['rates']['inf_engy'] = np.array([], dtype='float64')

    ######## -------- Formats Output -------- ########

    out = {}
    out['atomic'] = copy.deepcopy(FAC['atomic'])
    out['input_file'] = copy.deepcopy(FAC)
    out['rates'] = copy.deepcopy(FAC['rates'])

    # Output
    return out

############################################################
#
#               Transition Data Files
#
############################################################

# Reads energy level file
def _en(
    FAC=None,
    fil=None,
    ele=None,
    nele=None,
    Zele=None,
    ):

    # Useful constants
    eV2invcm = 8065.73 # [1/cm/eV]

    # Initializes output
    FAC['atomic'] = {}

    # Reads FAC energy levelfile
    en = rfac.read_en(
        fil+'a.en'
        )

    # Charge states
    FAC['atomic']['element'] = ele              # Species name
    FAC['atomic']['charge_state'] = Zele-nele   # Electronic charge
    FAC['atomic']['iz0'] = Zele                 # Nuclear charge
    FAC['atomic']['iz1'] = Zele-nele +1         # Spectroscopic charge

    # Initializes data arrays
    ion_pot = []                # [1/cm], dim(nion,), ionization potentials
    ion_term = []               # dim(nion,), ionization states
    ion_pot_lvl = []            # dim(nion,), ionization state index
    ion_L = []                  # dim(nion,), ionization state orbital angular momentum

    config = []                 # dim(ntran,), state configuration
    L = []                      # dim(ntran,), state L quantum number
    S = []                      # dim(ntran,), state S quantum number
    w = []                      # dim(ntran,), state J quantum number
    energy = []                 # dim(ntran,), state energy level wrt ground
    #zpla = []                   # dim(ntran,nion), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #zpla1 = []                  # dim(ntran,nion), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Convert FAC indices to form assumed by ColRadPy
    FAC['lvl_indices'] = {}
    FAC['lvl_indices']['FAC'] = np.array([], dtype='int')
    FAC['lvl_indices']['ColRadPy'] = np.array([], dtype='int')
    ind_e = 0
    ind_i = 0

    # Loop over blocks
    for blk in np.arange(len(en[1])):
        # If block of excited states
        if en[1][blk]['NELE'] == nele:
            # Loop over states
            for st in np.arange(len(en[1][blk]['ILEV'])):
                config.append(
                    en[1][blk]['sname'][st].decode("utf-8")
                    )

                L.append(
                    en[1][blk]['VNL'][st]%100
                    )

                w.append(
                    en[1][blk]['2J'][st]/2
                    )

                S.append(
                    2*abs(w[-1]-L[-1]) +1
                    )

                energy.append(
                    en[1][blk]['ENERGY'][st] * eV2invcm
                    )

                # Index conversion
                FAC['lvl_indices']['FAC'] = np.append(
                    FAC['lvl_indices']['FAC'],
                     en[1][blk]['ILEV'][st]
                     )
                ind_e += 1
                FAC['lvl_indices']['ColRadPy'] = np.append(
                    FAC['lvl_indices']['ColRadPy'], 
                    ind_e
                    )


        # If block of ionized states
        else:
            # Loop over states
            for st in np.arange(len(en[1][blk]['ILEV'])):
                ion_term.append(
                    en[1][blk]['sname'][st].decode("utf-8")
                    )

                ion_pot.append(
                    en[1][blk]['ENERGY'][st] * eV2invcm
                    )

                ion_L.append(
                    en[1][blk]['VNL'][st]%100
                    )

                # Index conversion
                FAC['lvl_indices']['FAC'] = np.append(
                    FAC['lvl_indices']['FAC'], 
                    en[1][blk]['ILEV'][st]
                    )
                ind_i += 1
                FAC['lvl_indices']['ColRadPy'] = np.append(
                    FAC['lvl_indices']['ColRadPy'], 
                    ind_i
                    )

                ion_pot_lvl.append(
                    ind_i
                    )


    # Stores data
    FAC['atomic']['ion_pot'] = np.asarray(ion_pot)
    FAC['atomic']['ion_term'] = np.asarray(ion_term)
    FAC['atomic']['config'] = np.asarray(config)
    FAC['atomic']['L'] = np.asarray(L)
    FAC['atomic']['S'] = np.asarray(S)
    FAC['atomic']['w'] = np.asarray(w)
    FAC['atomic']['energy'] = np.asarray(energy)
    #FAC['atomic']['zpla'] = np.asarray(zpla)
    #FAC['atomic']['zpla1'] = np.asarray(zpla1)
    FAC['atomic']['zpla'] = -1*np.ones((len(config), len(ion_pot)))
    FAC['atomic']['zpla1'] = -1*np.ones((len(config), len(ion_pot)))
    FAC['atomic']['ion_pot_lvl'] = np.asarray(ion_pot_lvl)
    FAC['atomic']['ion_L'] = np.asarray(ion_L)

    # Output
    return FAC

# Reads Einstein coefficients
def _tr(
    FAC=None,
    fil=None,
    nele=None,
    Zele=None,
    ):

    # Init output dictionary
    upr = []    # dim(ntrans,)
    lwr = []    # dim(ntrans,)
    a_val = []  # [1/s], dim(ntrans,)
    upr_ColRadPy = []
    lwr_ColRadPy = []

    # Reads transition rate data file
    tr = rfac.read_tr(fil+'a.tr')

    # Loop over blocks
    for blk in np.arange(len(tr[1])):
        # Loop over transitions
        for tran in np.arange(len(tr[1][blk]['lower_index'])):
            upr.append(
                tr[1][blk]['upper_index'][tran]
                )

            lwr.append(
                tr[1][blk]['lower_index'][tran]
                )

            a_val.append(
                tr[1][blk]['rate'][tran]
                )

            # Converts FAC indices to ColRadPy notation
            ind_upr = np.where(FAC['lvl_indices']['FAC'] == upr[-1])[0][0]
            ind_lwr = np.where(FAC['lvl_indices']['FAC'] == lwr[-1])[0][0]
            upr_ColRadPy.append(
                FAC['lvl_indices']['ColRadPy'][ind_upr]
                )
            lwr_ColRadPy.append(
                FAC['lvl_indices']['ColRadPy'][ind_lwr]
                )

    # Transition array from TRTable()
    trans_ColRadPy = np.vstack((upr_ColRadPy, lwr_ColRadPy)).T # dim(ntrans,2)

    # Includes two-photon emission
    #  Only works for H-like and He-like
    if nele <= 2:
        a_2photon = crm.TwoPhoton(Zele, nele-1) # [1/s]

        # 2s(0.5) for H-like
        if nele == 1:
            ind = np.where(FAC['atomic']['config'] == '2s1')[0][0]
        # 1s 2s (J=0) for He-like
        elif nele == 2:
            ind = np.where(
                (FAC['atomic']['config'] == '1s1.2s1')
                & (FAC['atomic']['w'] == 0)
                )[0][0]

        ind2 = np.where(
            np.all(
                trans_ColRadPy == np.asarray([
                    FAC['lvl_indices']['ColRadPy'][ind], 1
                    ]),
                axis = 1
                )
            )[0][0]

        a_val[ind2] += a_2photon

    ## NOTE: ColRadPy assumes the index of 'a_val' corresponds to 'col_transitions'
    # But, not all transitions have collisional exication or Einstein coefficients
    # Therefore need to fill
    col_trans = FAC['rates']['excit']['col_transitions'].copy()
    col_excit = FAC['rates']['excit']['col_excit'].copy()
    all_trans = np.unique(
        np.concatenate((
            trans_ColRadPy,
            col_trans,
            ), axis=0),
        axis=0)
    all_a_val = np.ones(all_trans.shape[0])*1e-30 # dim(ntrans,)
    all_excit = np.ones((all_trans.shape[0], col_excit.shape[1]))*1e-30 # dim(ntrans, ntemp)

    for tt in np.arange(all_trans.shape[0]):
        # Finds the indices for this transition
        inda = np.where(
            np.all(
                trans_ColRadPy == all_trans[tt,:],
                axis = 1
                )
            )[0]
        indc = np.where(
            np.all(
                col_trans == all_trans[tt,:],
                axis = 1
                )
            )[0]

        # Error check
        if len(inda) != 1 and len(indc) != 1:
            print(all_trans[tt,:])
            print('xxxx')

        # Fills arrays
        if len(inda) == 1:
            all_a_val[tt] = a_val[inda[0]]

        if len(indc) == 1:
            all_excit[tt,:] = col_excit[indc[0],:]

    # Output
    FAC['rates']['a_val'] = all_a_val                       # [1/s], dim(ntrans,)
    FAC['rates']['excit']['col_transitions'] = all_trans    # dim(ntrans,2)
    FAC['rates']['excit']['col_excit'] = all_excit          # [upsilon], dim(ntrans, ntemp)

    # Output
    return FAC

############################################################
#
#               Data Loading & Prep
#
############################################################

# Reads dielectronic recombination strength data files
## NOTE: Skips autoionization rates right now!!!
def _ai(
    nele=None,
    EEDF=None,
    Te=None,        # [eV], dim(ntemp,)
    FAC=None,
    fil=None,
    verbose=None,
    ):

    # Error check if H-like
    if nele == 1:
        return FAC

    # Init output data
    bnd = []
    fre = []
    bnd_ColRadPy = []
    fre_ColRadPy = []
    dE = [] # [eV]
    DC = [] # [eV *cm2]

    # Loop over blocks
    for blk in np.arange(len(ai[1])):
        # Loop over transitions
        for tran in np.arange(len(ai[1][blk]['Delta E'])):
            # Stores data
            dE.append(
                ai[1][blk]['Delta E'][trn]
                )

            DC.append(
                ai[1][blk]['DC strength'][trn]*1e-20
                )

            bnd.append(
                ai[1][blk]['bound_index'][trn]
                )
            fre.append(
                ai[1][blk]['free_index'][trn]
                )

            # Converts FAC indices to ColRadPy notation
            ind_bnd = np.where(FAC['lvl_indices']['FAC'] == bnd[-1])[0][0]
            ind_fre = np.where(FAC['lvl_indices']['FAC'] == fre[-1])[0][0]
            bnd_ColRadPy.append(
                FAC['lvl_indices']['ColRadPy'][ind_bnd]
                )
            fre_ColRadPy.append(
                FAC['lvl_indices']['ColRadPy'][ind_fre]
                )

    # Calculates rate coefficient assuming delta function energy-dependent cross-section
    ratec = convolve_EEDF(
        EEDF=EEDF,
        Te=Te,
        XS=DC,
        engyXS=dE,
        DC_flag=True,
        use_rel=True,
        ).T # [cm3/s],dim(ntrans,ntemp)

    # Transition array from AITable()
    trans_ColRadPy = np.vstack((fre_ColRadPy, bnd_ColRadPy)).T # dim(ntrans,2)

    # Combining together rad recomb and DC recomb 
    rad_trans = FAC['rates']['recomb']['recomb_transitions'].copy()
    rad_recomb = FAC['rates']['recomb']['recomb_excit'].copy()
    all_trans = np.unique(
        np.concatenate((
            trans_ColRadPy,
            rad_trans,
            ), axis=0),
        axis=0)
    all_recomb = np.ones((all_trans.shape[0], rad_recomb.shape[1]))*1e-30 # dim(ntrans, ntemp)

    for tt in np.arange(all_trans.shape[0]):
        # Finds the indices for this transition
        ind_dc = np.where(
            np.all(
                trans_ColRadPy == all_trans[tt,:],
                axis = 1
                )
            )[0]
        ind_rr = np.where(
            np.all(
                rad_trans == all_trans[tt,:],
                axis = 1
                )
            )[0]

        # Error check
        if len(ind_rr) != 1 and len(ind_dc) != 1:
            print(all_trans[tt,:])
            print('xxxx')

        # Fills arrays
        if len(ind_rr) == 1:
            all_recomb[tt,:] += rad_recomb[ind_rr[0],:]

        if len(ind_dc) == 1:
            all_recomb[tt,:] += ratec[ind_dc[0],:]

    # Output
    FAC['rates']['recomb']['recomb_transitions'] = all_trans    # dim(ntrans,2)
    FAC['rates']['recomb']['recomb_excit'] = all_recomb         # [cm3/s], dim(ntrans, ntemp)

    return FAC

# Reads cross-section data files
def _get_xs(
    EEDF=None,
    Te=None,       # [eV], dim(ntemp,)
    FAC=None,
    fil=None,
    vebose = None,
    reacts = None,
    ):

    # Error check
    if 'ce' not in reacts:
        print('NEED TO INCLUDE COLLISIONAL EXCITATION DATA IN MODELING!!!')
        sys.exit(1)

    # Fill blank
    if 'rr' not in reacts:
        FAC['rates']['recomb'] = {}
        FAC['rates']['recomb']['recomb_transitions'] = np.asarray([])
        FAC['rates']['recomb']['recomb_excit'] = np.asarray([])
    if 'ci' not in reacts:
        FAC['rates']['ioniz'] = {}
        FAC['rates']['ioniz']['ion_transitions'] = np.asarray([])
        FAC['rates']['ioniz']['ion_excit'] = np.asarray([])

    # Useful constants
    eV2K = 11604.5

    # Loop over reactions
    for react in reacts:
        # Skip
        if react not in ['ce', 'rr', 'ci']:
            continue

        # Loads date from ascii file
        XSdata = utils._read_ascii(
            fil = fil,
            react = react,
            )

        # Saves temperature grid data
        if 'temp_grid' not in FAC.keys():
            # NOTE: ColRadPy assumes Te is in Kelvin within the data file (per ADF04 standard)
            FAC['temp_grid'] = Te*eV2K # [K], dim(nt,)

        # Performs EEDF convolution to cross-sections
        (
            data,
            st_upr, st_lwr,
            ratec, XS, engy
            ) = utils._conv_ascii2colradpy(
                FAC = FAC,
                XSdata = XSdata,
                EEDF = EEDF,
                Te = Te,
                react = react,
                verbose = verbose,
                )

        # Dictionary labeling
        if react == 'ce':
            lbl = 'excit'
            lbl_trans = 'col_transitions'
            lbl_excit = 'col_excit'
        elif react == 'rr':
            lbl = 'recomb'
            lbl_trans = 'recomb_transitions'
            lbl_excit = 'recomb_excit'
        elif react == 'ci':
            lbl = 'ioniz'
            lbl_trans = 'ion_transitions'
            lbl_excit = 'ion_excit'

        # Formats output
        FAC['rates'][lbl] = {}

        # Stores transitions arrays
        if react in ['ce', 'rr']:
            # If excit -- (upr,lwr) states in ColRadPy indices
            # If recomb -- Z+1 state -> Z state
            FAC['rates'][lbl][lbl_trans] = np.vstack(
                (np.asarray(st_upr), np.asarray(st_lwr))
                ).T # dim(ntrans,2)
        else:
            # If ioniz -- Z state -> Z+1 state
            FAC['rates'][lbl][lbl_trans] = np.vstack(
                (np.asarray(st_lwr), np.asarray(st_upr))
                ).T # dim(ntrans,2)

        # Stores rate coefficient data
        # If excit -- [upsilon], if recomb -- [cm3/s], if ioniz -- [reduced]
        FAC['rates'][lbl][lbl_excit] = np.asarray(data) # dim(ntrans, nt)

        if verbose == 1:
            FAC['rates'][lbl]['XS_cm2'] = np.asarray(XS) # dim(ntrans, nE), [cm2]
            FAC['rates'][lbl]['engy_eV'] = np.asarray(engy) # dim(ntrans, nE), [eV]
            if react == 'ce':
                # Stores transition array before padding to match Einstein coefficient
                FAC['rates'][lbl]['col_trans_unfill'] = FAC['rates'][lbl][lbl_trans].copy()
            if react != 'rr':
                # Stores rate coefficients in SI units
                FAC['rates'][lbl]['ratec_cm3/s'] = np.asarray(ratec) # dim(ntrans, nt), [cm3/s]

    # Output
    return FAC

# Reads Maxwell-averaged data files
def _get_mr(
    FAC = None,
    fil = None,
    verbose = None,
    reacts = None,
    ):

    # Error check
    if 'ce' not in reacts:
        print('NEED TO INCLUDE COLLISIONAL EXCITATION DATA IN MODELING!!!')
        sys.exit(1)

    # Fill blank
    if 'rr' not in reacts:
        FAC['rates']['recomb'] = {}
        FAC['rates']['recomb']['recomb_transitions'] = np.asarray([])
        FAC['rates']['recomb']['recomb_excit'] = np.asarray([])
    if 'ci' not in reacts:
        FAC['rates']['ioniz'] = {}
        FAC['rates']['ioniz']['ion_transitions'] = np.asarray([])
        FAC['rates']['ioniz']['ion_excit'] = np.asarray([])

    # Useful constants
    eV2K = 11604.5

    # Loop over reactions
    for react in reacts:
        # Skip
        if react not in ['ce', 'rr', 'ci']:
            continue

        # Reads data file
        mr = utils._read_mr(
            fil = fil,
            react = react,
            )

        # Saves temperature grid data
        if 'temp_grid' not in FAC.keys():
            # NOTE: ColRadPy assumes Te is in Kelvin within the data file (per ADF04 standard)
            FAC['temp_grid'] = np.asarray(mr['Te_eV'])*eV2K # [K], dim(nt,)

        # Converts data file to ColRadPy form
        data, st_upr, st_lwr, ratec = utils._conv_mr2colradpy(
            FAC = FAC,
            mr = mr,
            react = react,
            verbose = verbose,
            )

        # Dictionary labeling
        if react == 'ce':
            lbl = 'excit'
            lbl_trans = 'col_transitions'
            lbl_excit = 'col_excit'
        elif react == 'rr':
            lbl = 'recomb'
            lbl_trans = 'recomb_transitions'
            lbl_excit = 'recomb_excit'
        elif react == 'ci':
            lbl = 'ioniz'
            lbl_trans = 'ion_transitions'
            lbl_excit = 'ion_excit'

        # Formats output
        FAC['rates'][lbl] = {}

        # Stores transitions arrays
        if react in ['ce', 'rr']:
            # If excit -- (upr,lwr) states in ColRadPy indices
            # If recomb -- Z+1 state -> Z state
            FAC['rates'][lbl][lbl_trans] = np.vstack(
                (np.asarray(st_upr), np.asarray(st_lwr))
                ).T # dim(ntrans,2)
        else:
            # If ioniz -- Z state -> Z+1 state
            FAC['rates'][lbl][lbl_trans] = np.vstack(
                (np.asarray(st_lwr), np.asarray(st_upr))
                ).T # dim(ntrans,2)

        # Stores rate coefficient data
        # If excit -- [upsilon], if recomb -- [cm3/s], if ioniz -- [reduced]
        FAC['rates'][lbl][lbl_excit] = np.asarray(data) # dim(ntrans, nt)

        if verbose == 1:
            if react == 'ce':
                # Stores transition array before padding to match Einstein coefficient
                FAC['rates'][lbl]['col_trans_unfill'] = FAC['rates'][lbl][lbl_trans].copy()
            if react != 'rr':
                # Stores rate coefficients in SI units
                FAC['rates'][lbl]['ratec_cm3/s'] = np.asarray(ratec) # dim(ntrans, nt), [cm3/s]

    # Output
    return FAC