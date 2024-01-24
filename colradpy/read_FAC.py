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
    5) 2-photon emission data

'''

# Module
from pfac import rfac, crm
import os
import numpy as np
import copy
from colradpy.convolve_EEDF import convolve_EEDF as cE

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
    physics = None,     # if None -> looks for all file suffixes
    Te = None,          # if not using MaxwellRate files
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

    # FAC data file suffixes to search for
    if physics == 'incl_all':
        # If already Maxwell-averaged rate coefficients
        if use_mr:
            physics = [
                'en',       # ASCII-format Energy levels
                'tr',       # ASCII-format Einstein coefficients
                'ce.mr',    # Collisional excitation
                'rr.mr',    # Radiative recombination
                #'ai',    # Autoionization/dielectronic recombination
                'ci.mr',    # Collision ionization
                ]
        # If using cross-section files
        else:
            physics = [
                'en',      # ASCII-format Energy levels
                'tr',      # ASCII-format Einstein coefficients
                'ce',      # ASCII-format Collisional excitation
                #'rr',      # ASCII-format Radiative recombination
                #'ai',     # ASCII-format Autoionization/dielectronic recombination
                #'ci',      # ASCII-format Collision ionization
                ]

    ######## -------- Reads data -------- ########

    # Initialize output
    FAC = {}
    FAC['rates'] = {}

    # Energy levels
    if 'en' in physics:
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

    # Collisional excitation
    if 'ce' in physics:
        FAC = _ce(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            Te=Te,
            )
    elif 'ce.mr' in physics:
        FAC = _ce_mr(
            FAC=FAC,
            fil=fil,
            )
    # Error check
    else:
        print('NEED TO INCLUDE COLLISIONAL EXCITATION DATA IN MODELING!!!')
        sys.exit(1)

    # Einstein coefficients
    if 'tr' in physics:
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

    # Radiative recombination
    if 'rr' in physics:
        FAC = _rr(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            )
    elif 'rr.mr' in physics:
        FAC = _rr_mr(
            FAC=FAC,
            fil=fil,
            )
    # If empty
    else:
        FAC['rates']['recomb'] = {}
        FAC['rates']['recomb']['recomb_transitions'] = np.asarray([])
        FAC['rates']['recomb']['recomb_excit'] = np.asarray([])

    # Autoionization/dielectronic recombination
    if 'ai' in physics:
        FAC = _ai(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            )
    elif 'ai.mr' in physics:
        FAC = _ai_mr(
            FAC=FAC,
            fil=fil,
            )

    # Collisional ionization
    if 'ci' in physics:
        FAC = _ai(
            FAC=FAC,
            fil=fil,
            EEDF=EEDF,
            )
    elif 'ci.mr' in physics:
        FAC = _ci_mr(
            FAC=FAC,
            fil=fil,
            )
    # If empty
    else:
        FAC['rates']['ioniz'] = {}
        FAC['rates']['ioniz']['ion_transitions'] = np.asarray([])
        FAC['rates']['ioniz']['ion_excit'] = np.asarray([])

    # Charge exchange
    if 'cx' in physics:
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
#                      File Reading
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

# Reads collisional excitation cross-section data files
def _ce(
    EEDF=None,
    Te=None,       # [eV], dim(ntemp,)
    FAC=None,
    fil=None,
    trans_FAC=None,
    vebose = 1,
    ):

    # Useful constants
    eV2K = 11604.5

    # Saves temperature grid data
    # NOTE: ColRadPy assumes Te is in Kelvin within the data file (per ADF04 standard)
    FAC['temp_grid'] = Te*eV2K # [K], dim(nt,)

    # Initializes rate data
    data = np.zeros((trans_FAC.shape[0], len(Te))) # dim(ntrans,nt)
    trans_ColRadPy = np.zeros(trans_FAC.shape, dtype='int')
    XS = None

    # Reads FAC data file
    ce = rfac.read_ce(fil+'a.ce')

    # Loop over blocks
    for blk in np.arange(len(ce[1])):
        # Loop over transitions in block
        for trn in np.arange(len(ce[1][blk]['Delta E'])):
            # Init data storage
            if XS is None:
                XS = np.zeros((len(ce[1][blk]['EGRID']),trans_FAC.shape[0])) # dim(nE, ntrans)
                engyXS = np.zeros((len(ce[1][blk]['EGRID']),trans_FAC.shape[0])) # dim(nE, ntrans)
                dE = np.zeros(trans_FAC.shape[0])
                Bethe = np.zeros((trans_FAC.shape[0], 2))
                w_upr = np.zeros(trans_FAC.shape[0])


            # Upper and lower level indices
            upr = ce[1][blk]['upper_index'][trn]
            lwr = ce[1][blk]['lower_index'][trn]

            # Converts indices
            ind_upr = np.where(FAC['lvl_indices']['FAC'] == upr)
            ind_lwr = np.where(FAC['lvl_indices']['FAC'] == lwr)

            try:
                ind_tt = np.where(
                    (trans_FAC[:,0] == upr)
                    & (trans_FAC[:,1] == lwr)
                    )[0][0]
            except:
                print(upr)
                print(lwr)
                print('xxxxx')
                continue

            #import pdb
            #pdb.set_trace()

            trans_ColRadPy[ind_tt,0] = FAC['lvl_indices']['ColRadPy'][ind_upr]
            trans_ColRadPy[ind_tt,1] = FAC['lvl_indices']['ColRadPy'][ind_lwr]

            # Stores data
            XS[:,ind_tt] = ce[1][blk]['crosssection'][trn,:]*1e-20
            engyXS[:,ind_tt] = ce[1][blk]['EGRID']
            dE[ind_tt] = ce[1][blk]['TE0']
            Bethe[ind_tt,0] = ce[1][blk]['bethe'][trn]
            Bethe[ind_tt,1] = ce[1][blk]['born'][trn,0]
            w_upr[ind_tt] = FAC['atomic']['w'][ind_upr]

    # Calculates Rate coefficient data, [cm3/s], dim(ntrans, ntemp)
    if verbose == 1:
        (
            data, 
            )= cE(
            EEDF = EEDF,
            Te=Te,
            XS = XS,
            engyXS = engyXS,
            m = ce[1][blk]['ETYPE'],
            dE = dE,
            Bethe = Bethe,
            #Bethe = None,
            w_upr = w_upr,
            verbose=verbose,
            )
        data = data.T

    else:
        data = cE(
            EEDF = EEDF,
            Te=Te,
            XS = XS,
            engyXS = engyXS,
            m = ce[1][blk]['ETYPE'],
            dE = dE,
            Bethe = Bethe,
            #Bethe = None,
            w_upr = w_upr,
            verbose=verbose,
            ).T

    #import pdb
    #pdb.set_trace()
    if verbose == 1:
        ratec = data.copy()


    # Converts data to expected reduced form (Upsilons per ADF04 standard)
    for tt in np.arange(trans_FAC.shape[0]):
        data[tt,:] = _conv_rate2upsilon(
            data = data[tt,:],
            Te_eV = Te,
            ind_upr = trans_ColRadPy[tt,0] -1,
            ind_lwr = trans_ColRadPy[tt,1] -1,
            FAC = FAC,
            )

            
    # Formats output
    FAC['rates']['excit'] = {}
    FAC['rates']['excit']['col_transitions'] = trans_ColRadPy    # dim(ntrans,2), (upr,lwr) states in ColRadPy indices
    FAC['rates']['excit']['col_excit'] = data           # dim(ntrans,nt), [upsilon]
    if verbose == 1:
        FAC['rates']['excit']['col_trans_unfill'] = FAC['rates']['excit']['col_transitions'].copy()
        FAC['rates']['excit']['ratec_cm3/s'] = ratec # dim(ntrans,nt), [cm3/s]
        FAC['rates']['excit']['Bethe'] = Bethe # dim(ntrans, 2)

    # Output 
    return FAC

# Reads Maxwell-averaged collisional excitation data files
def _ce_mr(
    FAC=None,
    fil=None,
    verbose = 1,
    ):

    # Useful constants
    eV2K = 11604.5

    # Reads data file
    mr = _read_mr(
        fil=fil,
        data='ce'
        )

    # Saves temperature grid data
    # NOTE: ColRadPy assumes Te is in Kelvin within the data file (per ADF04 standard)
    FAC['temp_grid'] = np.asarray(mr['Te_eV'])*eV2K # [K], dim(nt,)

    # Init output
    st_lwr = []
    st_upr = []
    ratec = []

    # Loop over lower states
    for lwr in mr.keys():
        # Skips temp grid
        if lwr == 'Te_eV':
            continue

        # Converts indices
        ind_lwr = np.where(FAC['lvl_indices']['FAC'] == lwr)[0][0]
        
        # Loop over upper states
        for upr in mr[lwr]['coll_excit'].keys():
            # Converts indices
            ind_upr = np.where(FAC['lvl_indices']['FAC'] == upr)[0][0]
            st_lwr.append(
                FAC['lvl_indices']['ColRadPy'][ind_lwr]
                )
            st_upr.append(
                FAC['lvl_indices']['ColRadPy'][ind_upr]
                )

            ratec.append(
                mr[lwr]['coll_excit'][upr]
                ) # [cm3/s]

    # Convert to Upsilon form
    data = np.zeros((len(ratec), len(mr['Te_eV']))) # dim(ntrans, nt)
    for tt in np.arange(len(ratec)):
        data[tt,:] = _conv_rate2upsilon(
            data = ratec[tt],
            Te_eV = mr['Te_eV'],
            ind_upr = st_upr[tt] -1,
            ind_lwr = st_lwr[tt] -1,
            FAC = FAC,
            )

    # Formats output
    FAC['rates']['excit'] = {}
    FAC['rates']['excit']['col_transitions'] = np.vstack(
        (np.asarray(st_upr), np.asarray(st_lwr))
        ).T # dim(ntrans,2), (upr,lwr) states in ColRadPy indices
    FAC['rates']['excit']['col_excit'] = data # dim(ntrans, nt), [cm3/s] 
    if verbose == 1:
        FAC['rates']['excit']['col_trans_unfill'] = FAC['rates']['excit']['col_transitions'].copy()
        FAC['rates']['excit']['ratec_cm3/s'] = np.asarray(ratec) # dim(ntrans, nt), [cm3/s]

    # Output
    return FAC

# Reads Maxwell-averaged Radiative recombination data
def _rr_mr(
    FAC=None,
    fil=None,
    ):

    # Reads data file
    mr = _read_mr(
        fil=fil,
        data='rr'
        )

    # Init output
    state = []
    ion = []
    data = []

    # Loop over states with charge Z
    for st in mr.keys():
        # Skips temp grid
        if st == 'Te_eV':
            continue

        # Converts indices
        ind_st = np.where(FAC['lvl_indices']['FAC'] == st)[0][0]

        # Loop over states with charge Z+1
        for ii in mr[st]['rad_recomb'].keys():
            # Converts indices
            ind_ion = np.where(FAC['lvl_indices']['FAC'] == ii)[0][0]
            state.append(
                FAC['lvl_indices']['ColRadPy'][ind_st]
                )
            ion.append(
                FAC['lvl_indices']['ColRadPy'][ind_ion]
                )

            data.append(
                mr[st]['rad_recomb'][ii]
                ) # [cm3/s]

    # Formats output
    FAC['rates']['recomb'] = {}
    FAC['rates']['recomb']['recomb_transitions'] = np.vstack(
        (np.asarray(ion), np.asarray(state))
        ).T # dim(ntrans,2), Z+1 state -> Z state
    FAC['rates']['recomb']['recomb_excit'] = np.asarray(
        data
        ) # dim(ntrans, nt), [cm3/s] 

    # Ouput
    return FAC
        
# Reads Maxwell-averaged collisional ionization data
def _ci_mr(
    FAC=None,
    fil=None,
    ):

    # Reads data file
    mr = _read_mr(
        fil=fil,
        data='ci'
        )

    # Init output
    state = []
    ion = []
    data = []

    # Loop over states with charge Z
    for st in mr.keys():
        # Skips temp grid
        if st == 'Te_eV':
            continue

        # Converts indices
        ind_st = np.where(FAC['lvl_indices']['FAC'] == st)[0][0]

        # Loop over states with charge Z+1
        for ii in mr[st]['coll_ion'].keys():
            # Converts indices
            ind_ion = np.where(FAC['lvl_indices']['FAC'] == ii)[0][0]
            state.append(
                FAC['lvl_indices']['ColRadPy'][ind_st]
                )
            ion.append(
                FAC['lvl_indices']['ColRadPy'][ind_ion]
                )

            data.append(
                mr[st]['coll_ion'][ii]
                ) # [cm3/s]

            # Converts ionization data to adf04 reduced form
            data[-1] = _conv_rate2reduct(
                data = data[-1],
                Te_eV = mr['Te_eV'],
                ind_st = state[-1] -1,
                ind_ion = ion[-1] -1,
                FAC = FAC,
                )

    # Formats output
    FAC['rates']['ioniz'] = {}
    FAC['rates']['ioniz']['ion_transitions'] = np.vstack(
        (np.asarray(state), np.asarray(ion))
        ).T # dim(ntrans,2), Z state -> Z+1 state
    FAC['rates']['ioniz']['ion_excit'] = np.asarray(
        data
        ) # dim(ntrans, nt), [reduced rate] 

    # Ouput
    return FAC

############################################################
#
#                      Utilities
#
############################################################

# Converts rate coefficeint to adf04 reduced form
def _conv_rate2reduct(
    data = None,    # [cm3/s], dim(ntemp,), list
    Te_eV = None,   # [eV], dim(ntemp,)
    ind_st = None,
    ind_ion = None,
    FAC = None,
    ):

    # Useful constants
    eV2invcm = 8065.73 # [1/cm/eV]

    # Init output
    out = []

    for tt in np.arange(len(Te_eV)):
        # Calculates reducted rate
        out.append(
            data[tt] * np.exp(
                (FAC['atomic']['ion_pot'][ind_ion] - FAC['atomic']['energy'][ind_st])
                / eV2invcm
                / Te_eV[tt]
                )
            )
    
    # Output
    return out

# Converts rate coefficient to adf04 Upsilon form
def _conv_rate2upsilon(
    data = None,    # [cm3/s], dim(ntemp,), array
    Te_eV = None,   # [eV], dim(ntemp,)
    ind_upr = None,
    ind_lwr = None,
    FAC = None,
    ):

    # Useful constants
    eV2invcm = 8065.73 # [1/cm/eV]

    return data * (
            np.sqrt(np.asarray(Te_eV)/13.6058)
            /2.1716e-8
            *(1 + 2*FAC['atomic']['w'][ind_upr])
            * np.exp(
                abs(
                    FAC['atomic']['energy'][ind_upr]
                    -FAC['atomic']['energy'][ind_lwr]
                    )
                /(np.asarray(Te_eV)*eV2invcm)
                )
            )

# Reads cross-section data files
def _read_ascii(
    fil = None,
    data = None,
    ):

    # Reads data file
    if data == 'ce':
        data_fil = rfac.read_ce(fil+'a.ce')
        lwr_lbl = 'lower_index'
        upr_lbl = 'upper_index'
        data_lbl = 'crosssection'
    elif data == 'rr':
        data_fil = rfac.read_rr(fil+'a.ce')
        lwr_lbl = 'bound_index'
        upr_lbl = 'free_index'
        data_lbl = 'RR crosssection'
    elif data == 'ci':
        data_fil = rfac.read_rr(fil+'a.ce')
        lwr_lbl = 'bound_index'
        upr_lbl = 'free_index'
        data_lbl = 'crosssection'

    # Initializes output dictionary
    out = {}

    # Loop over blocks
    for blk in np.arange(len(data_fil[1])):
        # Loop over transitions
        for trn in np.arange(len(data_fil[1][blk]['Delta E'])):
            # Transition indices
            lwr = data_fil[1][blk][lwr_lbl][trn]
            upr = data_fil[1][blk][upr_lbl][trn]

            # If new lower state
            if lwr not in out.keys():
                out[lwr] = {}
            out[lwr][upr] = {}

            # Stores cross-section data, [cm2]
            out[lwr][upr]['XS'] = data_fil[1][blk][data_lbl][trn,:]*1e-20

            # Stores energy grid in terms of incident electron energy, [eV]
            if data_fil[1][blk]['ETYPE'] == 0:
                out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID']
            elif data_fil[1][blk]['ETYPE'] == 1:
                out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID'] + data_fil[1][blk]['TE0']

            # Stores parameters for high-energy behavior
            if data == 'ce':
                out[lwr][upr]['limit'] = np.asarray([
                    data_fil[1][blk]['bethe'][trn],
                    data_fil[1][blk]['born'][trn,0]
                    ])

    # Output
    return out

# Read Maxwellian-averaged data files
def _read_mr(
    fil = None,
    data = None,
    ):

    # Reads data file
    f = open(
        fil+data+'.mr',
        'r'
        )

    # Initializes output dictionary
    out = {}

    # Data orginaization labels
    if data == 'ce':
        label = 'coll_excit'
    elif data == 'rr':
        label = 'rad_recomb'
    elif data == 'ci':
        label = 'coll_ion'

    # Loop over lines
    for line in f:
        # Skip line breaks
        if line == '\n':
            continue

        # If reading a header
        if line.split(' ')[0] == '#':
            # Lower level
            lwr = int(line.split('\t')[0][1:])
            # Upper level
            upr = int(line.split('\t')[2])

            # Number of temperature points
            ntemp = int(line.split('\t')[5])
            indt = 0

            # If temperature mesh not yet included
            if 'Te_eV' not in out.keys():
                out['Te_eV'] = []

            # If new lower state
            if lwr not in out.keys():
                out[lwr] = {}
                out[lwr][label] = {}

            out[lwr][label][upr] = []

        # If reading data
        else:
            line = line.replace('\t', ' ')
            # If need to add temperature mesh
            if len(out['Te_eV']) < ntemp:
                out['Te_eV'].append(
                    float(line.split(' ')[1])
                    )

            # Adds rate coefficient data, [cm3/s]
            out[lwr][label][upr].append(
                float(line.replace('  ', ' ').split(' ')[-1])*1e-10
                )

            # Increases index to next temp point
            indt += 1

    # Output
    return out
