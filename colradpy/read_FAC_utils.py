'''

read_FAC_utils.py is a script that houses general-purpose
functions for reading FAC data

cjperks
Jan 24, 2024


'''

# Module
from pfac import rfac
import os
import numpy as np
import copy
from colradpy.convolve_EEDF import convolve_EEDF

############################################################
#
#            Converting Data to ColRadPy format
#
############################################################

# Converts data read from _read_ascii to form wanted by ColRadPy
# NOTE: Performs EEDF convolution from cross-section to rate coeffs
def _conv_ascii2colradpy(
    FAC = None,
    XSdata = None,
    EEDF = None,
    Te = None,
    react = None,
    verbose = None,
    ):

    # Init output
    st_lwr = []
    st_upr = []
    data = []
    ratec = []
    XS = []
    engy = []

    # Loop over lower states
    for lwr in XSdata.keys():
        # Converts indices
        ind_lwr = np.where(FAC['lvl_indices']['FAC'] == lwr)[0][0]
        
        # Loop over upper states
        for upr in XSdata[lwr].keys():
            # Converts indices
            ind_upr = np.where(FAC['lvl_indices']['FAC'] == upr)[0][0]
            st_lwr.append(
                FAC['lvl_indices']['ColRadPy'][ind_lwr]
                )
            st_upr.append(
                FAC['lvl_indices']['ColRadPy'][ind_upr]
                )

            # Ionized state oribtal angular momentum for rr cross-section extrapolation
            if react == 'rr':
                ion_L = FAC['atomic']['ion_L'][st_upr[-1] -1]
            else:
                ion_L = 0

            # Calculates Rate coefficient data, [cm3/s], dim(ntrans, ntemp)
            if verbose == 1:
                (
                    tmp_data, 
                    tmp_XS,
                    tmp_engy
                    )= convolve_EEDF(
                    EEDF = EEDF,
                    Te=Te,
                    XS = XSdata[lwr][upr]['XS'],
                    engyXS = XSdata[lwr][upr]['engy'],
                    m = 0,
                    dE = np.asarray([XSdata[lwr][upr]['dE']]),
                    limit = XSdata[lwr][upr]['limit'][None,:],
                    w_upr = np.asarray([XSdata[lwr][upr]['w_upr']]),
                    w_lwr = np.asarray([XSdata[lwr][upr]['w_lwr']]),
                    ion_L = ion_L,
                    verbose=verbose,
                    use_rel = True,
                    react = react,
                    )
                data.append(tmp_data[:,0])
                XS.append(tmp_XS[:,0])
                engy.append(tmp_engy[:,0])

            else:
                data.append(
                    convolve_EEDF(
                        EEDF = EEDF,
                        Te=Te,
                        XS = XSdata[lwr][upr]['XS'],
                        engyXS = XSdata[lwr][upr]['engy'],
                        m = 0,
                        dE = np.asarray([XSdata[lwr][upr]['dE']]),
                        Bethe = XSdata[lwr][upr]['limit'][None,:],
                        w_upr = np.asarray([XSdata[lwr][upr]['w_upr']]),
                        w_lwr = np.asarray([XSdata[lwr][upr]['w_lwr']]),
                        ion_L = ion_L,
                        verbose=verbose,
                        use_rel = True,
                        react = react,
                        )[:,0]
                    )

            if verbose == 1:
                ratec.append(data[-1].copy())

            if react == 'ce':
                # Convert to Upsilon form
                data[-1] = _conv_rate2upsilon(
                    data = data[-1],
                    Te_eV = Te,
                    ind_lwr = st_lwr[-1] -1,
                    ind_upr = st_upr[-1] -1,
                    FAC = FAC,
                    )

            elif react == 'ci':
                # Converts ionization data to adf04 reduced form
                data[-1] = _conv_rate2reduct(
                    data = data[-1],
                    Te_eV = Te,
                    ind_st  = st_lwr[-1] -1,
                    ind_ion = st_upr[-1] -1,
                    FAC = FAC,
                    )

    # Output
    return data, st_upr, st_lwr, ratec, XS, engy

# Converts data read from _read_mr to form wanted by ColRadPy
def _conv_mr2colradpy(
    FAC = None,
    mr = None,
    react = None,
    verbose = None,
    ):

    # Init output
    st_lwr = []
    st_upr = []
    data = []
    ratec = []

    if react == 'ce':
        lbl = 'coll_excit'
    elif react == 'rr':
        lbl = 'rad_recomb'
    elif react == 'ci':
        lbl = 'coll_ion'

    # Loop over states with charge Z
    for lwr in mr.keys():
        # Skips temp grid
        if lwr == 'Te_eV':
            continue

        # Converts indices
        ind_lwr = np.where(FAC['lvl_indices']['FAC'] == lwr)[0][0]

        # Loop over states with charge Z+1
        for upr in mr[lwr][lbl].keys():
            # Converts indices
            ind_upr = np.where(FAC['lvl_indices']['FAC'] == upr)[0][0]
            st_lwr.append(
                FAC['lvl_indices']['ColRadPy'][ind_lwr]
                )
            st_upr.append(
                FAC['lvl_indices']['ColRadPy'][ind_upr]
                )

            data.append(
                mr[lwr][lbl][upr]
                ) # [cm3/s]

            if verbose == 1:
                ratec.append(data[-1].copy())

            if react == 'ce':
                # Convert to Upsilon form
                data[-1] = _conv_rate2upsilon(
                    data = data[-1],
                    Te_eV = mr['Te_eV'],
                    ind_lwr = st_lwr[-1] -1,
                    ind_upr = st_upr[-1] -1,
                    FAC = FAC,
                    )

            elif react == 'ci':
                # Converts ionization data to adf04 reduced form
                data[-1] = _conv_rate2reduct(
                    data = data[-1],
                    Te_eV = mr['Te_eV'],
                    ind_st  = st_lwr[-1] -1,
                    ind_ion = st_upr[-1] -1,
                    FAC = FAC,
                    )

    # Output
    return data, st_upr, st_lwr, ratec

############################################################
#
#            Reformating Rate Coefficient Data
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
    ind_lwr = None,
    ind_upr = None,
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

############################################################
#
#            Reading Data Files
#
############################################################

# Reads cross-section data files
def _read_ascii(
    fil = None,
    react = None,
    ):

    # Reads data file
    if react == 'ce':
        data_fil = rfac.read_ce(fil+'a.ce')
        lwr_lbl = 'lower_index'
        upr_lbl = 'upper_index'
        data_lbl = 'crosssection'
    elif react == 'rr':
        data_fil = rfac.read_rr(fil+'a.ce')
        lwr_lbl = 'bound_index'
        upr_lbl = 'free_index'
        data_lbl = 'RR crosssection'
    elif react == 'ci':
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

            # Stores transition energy, [eV]
            out[lwr][upr]['dE'] = data_fil[1][blk]['Delta E'][trn]

            if react == 'ce':
                # Stores energy grid in terms of incident electron energy, [eV]
                if data_fil[1][blk]['ETYPE'] == 0:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID']
                elif data_fil[1][blk]['ETYPE'] == 1:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID'] + data_fil[1][blk]['TE0']

                # Stores total angular momentum
                out[lwr][upr]['w_upr'] = data_fil[1][blk]['upper_2J'][trn]/2
                out[lwr][upr]['w_lwr'] = data_fil[1][blk]['lower_2J'][trn]/2

                # Stores parameters for high-energy behavior
                out[lwr][upr]['limit'] = np.asarray([
                    data_fil[1][blk]['bethe'][trn],
                    data_fil[1][blk]['born'][trn,0]
                    ])

            elif react == 'rr':
                # Stores energy grid in terms of photo-electron energy, [eV]
                if data_fil[1][blk]['ETYPE'] == 0:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID'] - data_fil[1][blk]['Delta E'][trn]
                elif data_fil[1][blk]['ETYPE'] == 1:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID'] 

                # Stores total angular momentum
                out[lwr][upr]['w_upr'] = data_fil[1][blk]['free_2J'][trn]/2
                out[lwr][upr]['w_lwr'] = data_fil[1][blk]['bound_2J'][trn]/2

                # Stores parameters for high-energy behavior, dim(4,)
                out[lwr][upr]['limit'] = data_fil[1][blk]['parameters'][trn,:]

            elif react == 'ci':
                # Stores energy grid in terms of incident electron energy, [eV]
                if data_fil[1][blk]['ETYPE'] == 0:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID']
                elif data_fil[1][blk]['ETYPE'] == 1:
                    out[lwr][upr]['engy'] = data_fil[1][blk]['EGRID'] + data_fil[1][blk]['Delta E'][trn]

                # Stores total angular momentum
                out[lwr][upr]['w_upr'] = data_fil[1][blk]['free_2J'][trn]/2
                out[lwr][upr]['w_lwr'] = data_fil[1][blk]['bound_2J'][trn]/2

                # Stores parameters for high-energy behavior, dim(4,)
                out[lwr][upr]['limit'] = data_fil[1][blk]['parameters'][trn,:]

    # Output
    return out

# Read Maxwellian-averaged data files
def _read_mr(
    fil = None,
    react = None,
    ):

    # Reads data file
    f = open(
        fil+data+'.mr',
        'r'
        )

    # Initializes output dictionary
    out = {}

    # Data orginaization labels
    if react == 'ce':
        label = 'coll_excit'
    elif react == 'rr':
        label = 'rad_recomb'
    elif react == 'ci':
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
