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
    fil = None,
    # Physics controls
    EEDF = None,        # if None -> assumes Maxwell-averages rates from pfac.fac.MaxwellRate
    physics = None,     # if None -> looks for all
    ):