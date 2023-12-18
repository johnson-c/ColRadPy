'''

read_FAC.py is a subfunction that reads atomic data created by
the Flexibly Atomic Code (FAC) to be used in collisional-
radiative modeling in ColRadPy

NOTE: Assumes you have a local version of FAC

cjperks
Dec 18, 2023

TO DO:
    1) Add non-Maxwellian convolution

'''

# Module
from pfac import rfac
import os
import numpy as np


#########################################
#
#           Main
#
#########################################

def read_FAC(
    # File Management
    fil = None,         # Common path to FAC files, excluding physics extentions
    # Physics controls
    EEDF = None,        # if None -> assumes Maxwell-averages rates from pfac.fac.MaxwellRate
    physics = None,     # if None -> looks for all file suffixes
    ):

    #### ---- Determines which files to search for ---- ####

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
                'a.en',     # ASCII-format Energy levels
                'a.tr',     # ASCII-format Einstein coefficients
                'ce.mr',    # Collisional excitation
                'rr.mr',    # Radiative recombination
                #'ai.mr'    # Autoionization/dielectronic recombination
                'ci.mr',    # Collision ionization
                ]
        # If using cross-section files
        else:
            physics = [
                'a.en',      # ASCII-format Energy levels
                'a.tr',      # ASCII-format Einstein coefficients
                'a.ce',      # ASCII-format Collisional excitation
                'a.rr',      # ASCII-format Radiative recombination
                #'a.ai'      # ASCII-format Autoionization/dielectronic recombination
                'a.ci',      # ASCII-format Collision ionization
                ]
