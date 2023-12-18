'''

read_FAC.py is a subfunction that reads atomic data created by
the Flexibly Atomic Code (FAC) to be used in collisional-
radiative modeling in ColRadPy

NOTE: Assumes you have a local version of FAC

cjperks
Dec 18, 2023

TO DO:
    1) Add non-Maxwellian convolution
    2) Missing energy level data

'''

# Module
from pfac import rfac
import os
import numpy as np


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
    ):

    ######## -------- Determines which files to search for -------- ########

    # Electron energy distribution settings
    if EEDF is None:
        # Use Maxwell-averaged rates from pfac.fac.MaxwellRate
        use_mr = True
    else:
        print('NON-MAXWELLIAN ELECTRONS NOT IMPLEMENTED YET')
        sys.exit(1)

    # FAC data file suffixes to search for
    if physics is None:
        # If already Maxwell-averaged rate coefficients
        if use_mr:
            physics = [
                'en',     # ASCII-format Energy levels
                'tr',     # ASCII-format Einstein coefficients
                'ce.mr',    # Collisional excitation
                'rr.mr',    # Radiative recombination
                #'ai.mr'    # Autoionization/dielectronic recombination
                'ci.mr',    # Collision ionization
                ]
        # If using cross-section files
        else:
            physics = [
                'en',      # ASCII-format Energy levels
                'tr',      # ASCII-format Einstein coefficients
                'ce',      # ASCII-format Collisional excitation
                'rr',      # ASCII-format Radiative recombination
                #'ai'      # ASCII-format Autoionization/dielectronic recombination
                'ci',      # ASCII-format Collision ionization
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

    # Einstein coefficients
    if 'tr' in physics:
        FAC, trans = _tr(
            FAC=FAC,
            fil=fil,
            )
    # Error check
    else:
        print('NEED TO INCLUDE EINSTEIN COEFFICIENT DATA IN MODELING!!!')
        sys.exit(1)

    # Collisional excitation
    if 'ce' in physics:
        FAC = _ce(
            FAC=FAC,
            fil=fil,
            trans=trans,
            EEDF=EEDF,
            )
    elif 'ce.mr' in physics:
        FAC = _ce_mr(
            FAC=FAC,
            fil=fil,
            trans=trans,
            )
    # Error check
    else:
        print('NEED TO INCLUDE COLLISIONAL EXCITATION DATA IN MODELING!!!')
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
            fil=fil
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
            fil=fil
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
            fil=fil
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

    ######## -------- Formats Output -------- ########

    # Output
    return FAC

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
    ev2invcm = 8065.73 # [/1cm/eV]

    # Initializes output
    FAC['atomic'] = {}

    # Reads FAC energy levelfile
    en = rfac.read_en(
        fil+'a.en'
        )

    # Charge states
    FAC['atomic']['element'] = ele          # Species name
    FAC['atomic']['charge_state'] = nele    # Electronic charge
    FAC['atomic']['iz0'] = Zele             # Nuclear charge
    FAC['atomic']['iz1'] = nele +1          # Spectroscopic charge

    # Initializes data arrays
    ion_pot = []                # [1/cm], dim(nion,), ionization potentials
    ion_term = []               # dim(nion,), ionization states
    ion_pot_lvl = []            # dim(nion,), ionization state index

    config = []                 # dim(ntran,), state configuration
    L = []                      # dim(ntran,), state L quantum number
    S = []                      # dim(ntran,), state S quantum number, !!!!!!!!!!!!!!!!
    w = []                      # dim(ntran,), state J quantum number
    energy = []                 # dim(ntran,), state energy level wrt ground
    zpla = []                   # dim(ntran,nion), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zpla1 = []                  # dim(ntran,nion), !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Loop over blocks
    for blk in en[1].keys():
        # If block of excited states
        if en[1][blk]['NELE'] == nele:
            # Loop over states
            for st in np.arange(len(en[1][blk]['ILEV'])):
                config.append(
                    en[1][blk]['sname'][st].decode("utf-8")
                    )

                L.append(
                    en[1][blk]['VNL'][st][-1]
                    )

                w.append(
                    en[1][blk]['2J'][st]/2
                    )

                energy.append(
                    en[1][blk]['ENERGY'][st] * eV2invcm
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

                ion_pot_lvl.append(
                    en[1][blk]['ILEV'][st] +1
                    )

    # Stores data
    FAC['atomic']['ion_pot'] = np.asarray(ion_pot)
    FAC['atomic']['ion_term'] = np.asarray(ion_term)
    FAC['atomic']['config'] = np.asarray(config)
    FAC['atomic']['L'] = np.asarray(L)
    FAC['atomic']['S'] = np.asarray(S)
    FAC['atomic']['w'] = np.asarray(w)
    FAC['atomic']['energy'] = np.asarray(energy)
    FAC['atomic']['zpla'] = np.asarray(zpla)
    FAC['atomic']['zpla1'] = np.asarray(zpla1)
    FAC['atomic']['ion_pot_lvl'] = np.asarray(ion_pot_lvl)

    # Output
    return FAC

# Reads Einstein coefficients
def _tr(
    FAC=None,
    fil=None,
    ):

    # Init output dictionary
    upr = []    # dim(ntrans,)
    lwr = []    # dim(ntrans,)
    a_val = []  # [1/s], dim(ntrans,)

    # Reads transition rate data file
    tr = rfac.read_tr(fil+'a.tr')

    # Loop over blocks
    for blk in np.arange(len(tr[1])):
        # Loop over transitions
        for tran in np.arange(len(tr[1][blk]['lower_index'])):
            upr.append(
                tr[1][blk]['upper_index'][tran] +1
                )

            lwr.append(
                tr[1][blk]['lower_index'][tran] +1
                )

            a_val.append(
                tr[1][blk]['rate'][tran]
                )

    # Formats output
    FAC['rates']['a_val'] = np.asarray(a_val)   # [1/s], dim(ntrans,)
    trans = np.vstack((upr,lwr)).T # dim(ntrans,2) -> form for coll. excit transitions

    # Output
    return FAC, trans

# Reads Maxwell-averaged collisional excitation data files
def _ce_mr(
    FAC=None,
    fil=None,
    trans=None,
    ):

    # Reads data file
    mr = _read_mr(
        fil=fil,
        data='ce'
        )

    # Saves temperature grid data
    FAC['temp_grid'] = np.asarray(mr['Te_eV']) # [eV], dim(nt,)

    # Initializes rate data
    data = np.zeros((trans.shape[0], len(mr['Te_eV']))) # dim(ntrans,nt)

    # Loop over transitions
    for tt in np.arange(trans.shape[0]):
        # Upper and lower level indices
        upr = int(trans[tt,0] -1)
        lwr = int(trans[tt,1] -1)

        # Saves transition rate coefficients, [cm3/s]
        data[tt,:] = np.asarray(
            data[lwr]['coll_excit'][upr]
            )

    # Formats output
    FAC['rates']['excit'] = {}
    FAC['rates']['excit']['col_transitions'] = trans    # dim(ntrans,2), (upr,lwr) states
    FAC['rates']['excit']['col_excit'] = data           # dim(ntrans,nt), [cm3/s]

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

        state.append(
            st+1
            )

        # Loop over states with charge Z+1
        for ii in mr[st]['rad_recomb'].keys():
            ion.append(
                ii+1
                )

            data.append(
                mr[st]['rad_recomb'][ii]
                ) # [cm3/s]

    # Formats output
    FAC['rates']['recomb'] = {}
    FAC['rates']['recomb']['recomb_transitions'] = np.vstack(
        (np.asarray(ion), np.asarray(st))
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

        state.append(
            st+1
            )

        # Loop over states with charge Z+1
        for ii in mr[st]['coll_ion'].keys():
            ion.append(
                ii+1
                )

            data.append(
                mr[st]['coll_ion'][ii]
                ) # [cm3/s]

    # Formats output
    FAC['rates']['ioniz'] = {}
    FAC['rates']['ioniz']['ion_transitions'] = np.vstack(
        (np.asarray(st), np.asarray(ion))
        ).T # dim(ntrans,2), Z state -> Z+1 state
    FAC['rates']['ioniz']['ion_excit'] = np.asarray(
        data
        ) # dim(ntrans, nt), [cm3/s] 

    # Ouput
    return FAC

############################################################
#
#                      Utilities
#
############################################################


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
                    float(line.split('  ')[0])
                    )

            # Adds rate coefficient data, [cm3/s]
            out[lwr][label][upr].append(
                float(line.split('  ')[-1])*1e-10
                )

            # Increases index to next temp point
            indt += 1

    # Output
    return out
