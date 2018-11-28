import numpy as np
from r8yip import *

def r8necip(IZ, ion_pot, energy, zpla,temperature_grid):
    """THis function calculates ECIP rates, adapted from ADAS code.

    :param IZ: charge state
    :type IZ: int

    :param ion_pot: ionization potential
    :type ion_pot: float

    :param energy: Energy levels
    :type energy: array

    :param zpla: ionization potential
    :type zpla: array

    :param temp_grid: Temperature grid to return ECIP rates on
    :type temp_grid: array

    :returns: array --ECIP rates on the temperature grid

    """
    
    ZETA = zpla
    CL = 2.3
    DXIPOW = 1.5
    TK2ATE = 1.5789e5
    R2GAM  = 2.17161e-8 
    D150   = 1.50e2
    CR2GAM = CL*R2GAM
    CALF = 3.30048e-24
    XI_arr = ion_pot - energy
    X = np.array([ 0.26356 , 1.41340 , 3.59643 , 7.08581 ,
                   12.64080])
    W = np.array([ 0.521756   , 0.398667   , 0.0759424 ,
                   0.00361176 , 0.00002337])
    r8necip = np.zeros((len(energy),len(temperature_grid)))
    for z in range(0, len(XI_arr)):
        for zz in range(0,len(temperature_grid)):
            TE = temperature_grid[zz]
            XI = XI_arr[z]
            Z     = float( IZ+1 )
            ATE   = TK2ATE/TE
            EN    = Z/np.sqrt(XI)
            Y     = XI*ATE
            IZETA = ZETA[z]
            J1    = 5
            AI    = 0.0

            for  J in range (0,J1):
                V     =  X[J]
                B     =  V/Y
                B1    =  B+1.0
                C     =  np.sqrt(B1)
                R     =  (1.25*EN*EN+0.25)/Z
                DELTA =  (Z/EN)*(R+2.0*EN*EN*C/((B+2.0)*Z*Z))/(C+np.sqrt(B))
                C1    =  1.0/(B+2.0)
                F     =  1.0
                C4    =  R8YIP(0.0,DELTA)
                C2    =  (C1*(B-B1*C1*np.log(B1))+0.65343*(1.0-1.0/(B1**3))*C4/EN)*F
                AI    =  AI+W[J]*C2
                Q     =  4.0*C2/B1
                
            ALFRED = CALF * (ATE ** 1.5) * 8.68811e-8*np.sqrt(ATE)*ZETA[z]*AI/XI

            if( Y < 150):
                R8NECIP = 8.68811e-8*np.sqrt(ATE)*np.exp(-Y)*ZETA*AI/XI
            else:
                R8NECIP = 0.0

            r8necip[z,zz] = ALFRED
    return r8necip
