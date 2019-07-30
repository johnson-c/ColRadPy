################################################################################
# file name         : r8yip.py
# author            : Curt Johnson
# description       : This code uses the ECIP forumula to calculate ionization
# version           : 0.1
# python version    : 2.7.12 ipython 2.4.1
# dependencies      : numpy
#
# This code was entirely based on the ADAS fortran routine r8ecip.for and r8yip.for.
# Origal fortran authors Jonathan Nash, Hugh Summers, Martin O'Mulane
# It has just been transplated to python. Variable names may have been changed.
# See the Summers appelton review for the theory
#
################################################################################


import numpy as np

def R8YIP(XI,DELTA):
    """This code was entirely based on the ADAS fortran routine r8ecip.for and r8yip.for.
       It has just been transplated to python. Variable names may have been changed.
       See the Summers appelton review for the theory

    :param XI: -
    :type XI: -

    :param DELTA: -
    :type DELTA:  - 

    """
    P1 = 5.0E0 / 30.0E0 
    A = np.array([ 2.3916E0  , 1.6742E0  , 1.2583E0  , 9.738E-1  ,
                        7.656E-1  , 6.078E-1  , 4.8561E-1 , 3.8976E-1 ,
                        3.1388E-1 , 2.5342E-1 , 2.0501E-1 , 1.6610E-1 ,
                        1.3476E-1 , 1.0944E-1 , 8.896E-2  , 7.237E-2  ,
                        5.8903E-2 , 4.7971E-2 , 3.9086E-2 , 3.1860E-2])
    B = np.array([ 1.0091E0  , 3.015E-1  , 1.314E-1  , 7.63E-2   ,
                        5.04E-2   , 3.561E-2  , 2.634E-2  , 1.997E-2  ,
                        1.542E-2  , 1.205E-2  , 9.50E-3   , 7.57E-3   ,
                        6.02E-3   , 4.84E-3   , 3.89E-3   , 3.12E-3   ,
                        2.535E-3  , 2.047E-3  , 1.659E-3  , 1.344E-3]) 
    W1 = 0.0E0
    W2 = 0.0E0
    W3 = 0.0E0

    XM = np.abs( XI )
    X  = XM + DELTA

    if (X <= 0.1E0):
        T = np.log( 1.1229E0 / X )
        Y = T + 0.25E0 * X**2 * ( 1.0E0 - 2.0E0 * T**2 )
    elif (X >= 2.0E0):
        T = 1.0E0 / X
        Y = 1.5707963E0 * ( 1.0E0 + T * ( 0.25E0 + T
                * ( -0.09375E0 + T * ( 0.0703125E0 - T * 0.0472E0 ) ) ) )
        W1 = -2.0E0 * X
    else:
        J  = int( 10.0E0 * X )
        X0 = 0.1E0 * float( J )
        R  = 10.0E0 * ( X - X0 )
        S  = 1.0E0 - R
        J=J-1
        Y  = R * A[J+1] + S * A[J] + P1 * ( R * ( R**2 - 1.0E0 )
                 * B[J+1] + S * ( S**2 - 1.0E0 ) * B[J] )
    if ((X <= 0.1E0) or  (X >= 2.0E0) or (XM > 0.01E0)):
        X  = XM**2 / ( XM + DELTA )
        X2 = 0.125E0 * XM**2 / ( XM + 3.0E0 * DELTA )
        R8YIP = Y / ( 1.0E0 + X2**2 )
        W2 = 3.1416E0 * XM - 0.9442E0 * X * ( 1.0E0 + X * ( 2.14E0 + 14.2E0 * X ) )/ ( 1.0E0 + X * ( 2.44E0 + 12.8E0 * X ) )
    else:
        R8YIP = Y
    if (XI < 0.0E0):
        W3 = 6.283E0 * XI

    W = W1 + W2 + W3
    R8YIP = R8YIP * np.exp( W )
    return R8YIP

def r8ecip(IZC, ion_pot, energy, zpla, temperature_grid):

    MXT = 5
    P1 = 157890.0
    X = np.array([ 0.26356E0  , 1.41340E0  , 3.59643E0  ,
                        7.08581E0  , 12.6408E0])
    W = np.array([ 5.21756E-1 , 3.98667E-1 , 7.59424E-2 ,
                        3.61176E-3 , 2.337E-5])
    Z   = float( IZC + 1 )    
    XI_arr = ion_pot - energy
    ecip = np.zeros((len(energy),len(temperature_grid)))
    for z in range(0,len(XI_arr)):
        for zz in range(0,len(temperature_grid)):
            TE = temperature_grid[zz]
            XI = XI_arr[z]* 9.11269E-06
            ZETA = zpla[z]
            ATE = P1 / TE
            EN  = Z / np.sqrt( XI )
            Y   = XI * ATE
            AI = 0.0E0
            
            for J in range(0, MXT-1):
                B  = X[J] / Y
                B1 = B + 1.0E0
                C  = np.sqrt( B1 )
                R  = ( 1.25E0 * EN**2 + 0.25E0 ) / Z

                D  = ( Z / EN ) * ( R + 2.0E0 * EN**2 * C / ( ( B + 2.0E0 ) * Z**2 ) ) / ( C + np.sqrt( B ) )

                C1 = 1.0E0 / ( B + 2.0E0 )
                F  = 1.0E0
                C4 = R8YIP( 0.0E0 , D )
                C2 = ( C1 * ( B - B1 * C1 * np.log( B1 ) ) + 0.65343E0 * ( 1.0E0 - 1.0E0 / B1**3 ) * C4 / EN ) * F
                AI = AI + W[J] * C2

            if (Y > 180.0E0):
                R8ECIP = 0.0E0
            else:
                R8ECIP = 8.68811E-8 * np.sqrt( ATE ) * np.exp( -Y ) * ZETA * AI / XI
            ecip[z,zz] = R8ECIP
    return ecip
