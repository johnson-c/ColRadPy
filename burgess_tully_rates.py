import numpy as np


def type1_xconvert(tconvert_grid,energy,c=1.5):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :returns: array --type1 burgess x values

    """
    return 1- np.log(c) / np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*\
                                               0.69488623/ (energy[0,:] - energy[1,:] ) + c))

def type1_yconvert(tconvert_grid,energy,col_excit,direct='f'):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param col_excit: excitation rates of the levels that are being transformed
    :type temp_grid: array

    :returns: array --type1 burgess y values

    """
    
    if(direct == 'b' or direct == 'B' or direct == 'back'):
        return col_excit * np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*\
                                    0.69488623/ (energy[0,:] - energy[1,:])  + np.exp(1)))
    else:
        return col_excit / np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*\
                                    0.69488623/ (energy[0,:] - energy[1,:])  + np.exp(1)))

def type2_xconvert(tconvert_grid,energy,c=1.5):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param c: scaling constant doesn't impact results
    :type c: array
    
    :returns: array --type2 burgess x values

    """
    
    return np.transpose(tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/\
                        (energy[0,:] - energy[1,:] )/( tconvert_grid.reshape(
                            len(tconvert_grid),1)*0.69488623/ (energy[0,:] - energy[1,:] ) +c))

def type2_yconvert(col_excit,direct='f'):

    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully

    :param col_excit: excitation rates of the levels that are being transformed
    :type temp_grid: array

    :returns: array --type2 burgess y values

    """


    
    return col_excit

def type3_xconvert(tconvert_grid,energy,c=1.5):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param c: scaling constant doesn't impact results
    :type c: array
    
    :returns: array --type3 burgess x values

    """
    
    return np.transpose((tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/ (energy[0,:] - energy[1,:] ))/\
        ( tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/ (energy[0,:] - energy[1,:]) +c))

def type3_yconvert(tconvert_grid,energy,col_excit,direct='f'):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param col_excit: excitation rates of the levels that are being transformed
    :type temp_grid: array

    :returns: array --type3 burgess y values

    """
    if(direct == 'b' or direct == 'B' or direct == 'back'):

        return col_excit/np.transpose(tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/\
                        (energy[0,:] - energy[1,:]) +1)
        
    else:
        return np.transpose(tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/\
                        (energy[0,:] - energy[1,:]) +1)*col_excit

def type4_xconvert(tconvert_grid,energy,c=1.5):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param c: scaling constant doesn't impact results
    :type c: array
    
    :returns: array --type4 burgess x values

    """
    
    return 1- np.log(c) / np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*0.69488623/\
                                               (energy[0,:] - energy[1,:]) + c))

def type4_yconvert(tconvert_grid,energy,col_excit,c=1.5,direct='f'):
    """This function extrapolates excitation rates through the burgess-tully method
       see On the anlysis of collision strengths and rate coefficients A burgess and J.A. Tully


    :param tconvert_grid: Temperature grid to do transform on
    :type energy: array


    :param energy: energies of the levels that are being transformed
    :type temp_grid: array

    :param col_excit: excitation rates of the levels that are being transformed
    :type temp_grid: array

    :returns: array --type4 burgess y values

    """
    if(direct == 'b' or direct == 'B' or direct == 'back'):
        return col_excit * np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*\
                                            0.69488623/ (energy[0,:] - energy[1,:])  + c))

    else:
        return col_excit / np.transpose(np.log( tconvert_grid.reshape(len(tconvert_grid),1)*\
                                            0.69488623/ (energy[0,:] - energy[1,:])  + c))


def lin_yconvert(coeff_lin_b,coeff_lin_m,xval_extrap):
    """ Some rates may not have infinite energy points in which case a
        linear extrapolation is done in the burgess tully space from the 
    last two calculated points

    :param coeff_lin_b: coefficients of 'b' in the y=mx+b equation
    :type energy: array

    :param coeff_lin_m: coefficients of 'm' in the y=mx+b equation
    :type temp_grid: array

    :param xval_extrap: xvalues for the extrapolated temperatures
    :type temp_grid: array

    :returns: array --linear extrapolation rate  values

    """    
    return np.transpose(np.exp( coeff_lin_b + coeff_lin_m*np.transpose(xval_extrap)))



def calc_coeffs_b(xval_arr,yval_arr,yval_inf):
    """ For an infinite energy extrapolation temperatures and excitation rates
    are converted to burgess-tully space. A linear extrapolation in the 
    brugess-tully space is done between the last calculated point and the 
    infinite energy point. This calculates the 'b' in the equation y = ax+b

    :param xval_arr: array of x values from calcualted points must have atleast 2
    :type energy: array

    :param yval_arr: array of y values from calcualted points must have atleast 2
    :type energy: array

    :param yval_inf: array of infinite energy y values from calcualtions
    :type energy: array

    :returns: array -- b coefficent
    """

    
    return (np.log(yval_inf) - np.log(yval_arr[:,-1])) / (1 -xval_arr[:,-1])#(np.log(yval_arr[:,-1]) -np.log(yval_inf))/(xval_arr[:,-1] -1)


def calc_coeffs_a(coeffs_b,xval_arr,yval_arr):
    """ For an infinite energy extrapolation temperatures and excitation rates
    are converted to burgess-tully space. A linear extrapolation in the 
    brugess-tully space is done between the last calculated point and the 
    infinite energy point. This calculates the 'a' in the equation y = ax+b

    :param xval_arr: array of x values from calcualted points must have atleast 2
    :type energy: array

    :param yval_arr: array of y values from calcualted points must have atleast 2
    :type energy: array

    :param yval_inf: array of infinite energy y values from calcualtions
    :type energy: array

    :returns: array -- a coefficent
    """    
    
    return np.log(yval_arr[:,-1])-coeffs_b*xval_arr[:,-1]#np.exp(np.log(yval_inf) -1*coeffs_b)

def calc_coeffs_lin_m(xval_arr,yval_arr):
    """ Some rates may not have infinite energy points in which case a
        linear extrapolation is done in the burgess tully space from the 
    last two calculated points this calculates the 'b' coefficient in the 
    equation y = mx+b

    :param xval_arr: array of x values from calcualted points must have atleast 2
    :type energy: array

    :param yval_arr: array of y values from calcualted points must have atleast 2
    :type energy: array

    :param yval_inf: array of infinite energy y values from calcualtions
    :type energy: array


    :returns: array -- m coefficent

    """

    
    return (np.log(yval_arr[:,np.shape(yval_arr)[1]-1]) -np.log(yval_arr[:,-2]))/\
        (xval_arr[:,np.shape(xval_arr)[1]-1] - xval_arr[:,-2])

def calc_coeffs_lin_b(xval_arr,yval_arr,coeff_lin_m):
    """ Some rates may not have infinite energy points in which case a
        linear extrapolation is done in the burgess tully space from the 
    last two calculated points this calculates the 'b' coefficient in the 
    equation y = mx+b

    :param xval_arr: array of x values from calcualted points must have atleast 2
    :type energy: array

    :param yval_arr: array of y values from calcualted points must have atleast 2
    :type energy: array

    :param yval_inf: array of infinite energy y values from calcualtions
    :type energy: array


    :returns: array -- b coefficent

    """
    
    return  np.log(yval_arr[:,-1]) - coeff_lin_m*xval_arr[:,-1]




def burgess_tully_rates(user_temp_grid, calc_temp_grid, col_transitions,col_excit,energy,w,a_val,S,L,inf_engy,c=1.5):

    """ Calculates excitation rates extrapolated to temperatures beyond the calculated
    rates in the file. If an infinite energy point in the file is present the extrapolation
    will be done between the last calculated point and the infinite energy point. A linear
    line is fit in the burgess tully space and extrapolated from there before being transformed
    back into temperature. If no infinite energy point is provided then an extrapolation is 
    provided from the last two calculated points by fitting a linear line in burgess tully
    space.

    :user_temp_grid: array of user defined temperatures
    :type user_temp_grid: array

    :calc_temp_grid: array of temperatures which were calculated in atomic code
    :type calc_temp_grid: array

    :col_transitions: 2d array of transitions (upper,lower)
    :type col_transitions: array

    :col_excit: 2d array of excitations between levels
    :type col_excit: array

    :energy: array of energies for the levels
    :type energy: array

    :w: array of w (j) for the levels
    :type w: array

    :a_val: array of einstien a coefficients
    :type a_val: array

    :a_val: array of einstien a coefficients
    :type a_val: array

    :S: array of spin values for the levels
    :type S: array

    :L: array of spin values for the levels
    :type L: array

    :inf_engy: array of infinite enery points for the transitions
    :type inf_engy: array

    :c: scale value for the extrapolation, just impacts how zoomed the plot is not actual results
    :type c: float

    :returns: array -- b coefficent

    """    



    burgtully_dict = {}
    burgtully_dict['burg_tully'] = {}
    burgtully_dict['burg_tully']['interp_temp_inds'] = np.where( user_temp_grid < np.max(calc_temp_grid)/11604.5)[0]
    burgtully_dict['burg_tully']['extrap_temp_inds'] = np.where( user_temp_grid > np.max(calc_temp_grid)/11604.5)[0]
    FBIG = 0.01
    FZERO = 1E-4
    ELU = np.abs(energy[col_transitions[:,0]-1] - energy[col_transitions[:,1]-1])/109737.26
    WTU = 2*w[col_transitions[:,0]-1]+1
    WTL = 2*w[col_transitions[:,1]-1]+1
    SS = 3.73491E-10*a_val*WTU/ELU**3
    FIN = 1/3.*ELU*SS/WTL
####################################################################################################
    #populate the indice arr. this lets us see which transistions correspond to which types
    #this is a list of array of the four types, the values are the indices of the transition
    burgtully_dict['burg_tully']['ind_arrs'] = []
    burgtully_dict['burg_tully']['ind_arrs'].append(np.where((FIN>FBIG) &(S[col_transitions[:,0]-1] ==
                                        S[col_transitions[:,1]-1]) & (np.abs(L[col_transitions[:,0]-1] -
                                                                    L[col_transitions[:,1]-1]) <=1))[0] )
    burgtully_dict['burg_tully']['ind_arrs'].append( np.where((S[col_transitions[:,0]-1] == \
                                                               S[col_transitions[:,1]-1]) &
                                                    ((np.abs(L[col_transitions[:,0]-1] -
                                                    L[col_transitions[:,1]-1]) >1) | (FIN<0.01)))[0] )

    burgtully_dict['burg_tully']['ind_arrs'].append(np.where(((FIN>0.01) | (FIN<FZERO)) &
                                        (S[col_transitions[:,0]-1] !=S[col_transitions[:,1]-1]) )[0])
                                         
    burgtully_dict['burg_tully']['ind_arrs'].append( np.where((FIN>FZERO) & (S[col_transitions[:,0]-1] !=
                                        S[col_transitions[:,1]-1]) & (FIN<FBIG))[0] )
    #this is the constant that sets the scale doesn't impact the actual results more for plotting
    burgtully_dict['burg_tully']['c'] = c

    #set up the x value and y value arrays for the transitions on the file temperature grid
    #these values are in burgess tully space
    burgtully_dict['burg_tully']['xval_arrs'] = []
    burgtully_dict['burg_tully']['yval_arrs'] = []

    type1_energies =  np.array([energy[col_transitions
                                [burgtully_dict['burg_tully']['ind_arrs'][0]][:,0] -1],energy[col_transitions[
                                    burgtully_dict['burg_tully']['ind_arrs'][0]][:,1] -1]])
    burgtully_dict['burg_tully']['xval_arrs'].append(type1_xconvert(calc_temp_grid, type1_energies,c))
    
    burgtully_dict['burg_tully']['yval_arrs'].append(type1_yconvert(calc_temp_grid,type1_energies,
                                        col_excit[burgtully_dict['burg_tully']['ind_arrs'][0]] ))

    type2_energies =  np.array([energy[col_transitions
                                [burgtully_dict['burg_tully']['ind_arrs'][1]][:,0] -1],energy[col_transitions[
                                    burgtully_dict['burg_tully']['ind_arrs'][1]][:,1] -1]])

    burgtully_dict['burg_tully']['xval_arrs'].append(type2_xconvert(calc_temp_grid, type2_energies,c))
    burgtully_dict['burg_tully']['yval_arrs'].append(type2_yconvert(\
                                                    col_excit[burgtully_dict['burg_tully']['ind_arrs'][1]]))


    type3_energies =  np.array([energy[col_transitions
                                [burgtully_dict['burg_tully']['ind_arrs'][2]][:,0] -1],energy[col_transitions[
                                    burgtully_dict['burg_tully']['ind_arrs'][2]][:,1] -1]])

    burgtully_dict['burg_tully']['xval_arrs'].append(type3_xconvert(calc_temp_grid,type3_energies,c))
    
    burgtully_dict['burg_tully']['yval_arrs'].append(type3_yconvert(calc_temp_grid,type3_energies,
                                                       col_excit[burgtully_dict['burg_tully']['ind_arrs'][2]]))


    type4_energies =  np.array([energy[col_transitions
                                [burgtully_dict['burg_tully']['ind_arrs'][3]][:,0] -1],energy[col_transitions[
                                    burgtully_dict['burg_tully']['ind_arrs'][3]][:,1] -1]])


    burgtully_dict['burg_tully']['xval_arrs'].append(type4_xconvert(calc_temp_grid,type4_energies,c))
    burgtully_dict['burg_tully']['yval_arrs'].append(type4_yconvert(calc_temp_grid,type4_energies,
                                                col_excit[burgtully_dict['burg_tully']['ind_arrs'][3]],c))

    #popoulated the infite energy points for the four types
    burgtully_dict['burg_tully']['yvals_inf'] = []
    burgtully_dict['burg_tully']['yvals_inf'].append( np.abs(\
                                            inf_engy[burgtully_dict['burg_tully']['ind_arrs'][0]]))
    burgtully_dict['burg_tully']['yvals_inf'].append( np.abs(\
                                            inf_engy[burgtully_dict['burg_tully']['ind_arrs'][1]]))
    burgtully_dict['burg_tully']['yvals_inf'].append( np.abs(\
                                            inf_engy[burgtully_dict['burg_tully']['ind_arrs'][2]]))
    burgtully_dict['burg_tully']['yvals_inf'].append( np.abs(\
                                            inf_engy[burgtully_dict['burg_tully']['ind_arrs'][3]]))

    burgtully_dict['burg_tully']['coeffs_b'] = []
    burgtully_dict['burg_tully']['coeffs_a'] = []
    burgtully_dict['burg_tully']['zero_inds'] = []
    burgtully_dict['burg_tully']['coeffs_lin_b'] = []
    burgtully_dict['burg_tully']['coeffs_lin_m'] = []    
####################################################################################################
    #type1 coeffseicients
    burgtully_dict['burg_tully']['coeffs_b'].append(calc_coeffs_b(burgtully_dict['burg_tully']['xval_arrs'][0],
                                                    burgtully_dict['burg_tully']['yval_arrs'][0],
                                                    burgtully_dict['burg_tully']['yvals_inf'][0]))

    burgtully_dict['burg_tully']['coeffs_a'].append(calc_coeffs_a(burgtully_dict['burg_tully']['coeffs_b'][0],
                                                                  burgtully_dict['burg_tully']['xval_arrs'][0],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][0]))

    burgtully_dict['burg_tully']['zero_inds'].append( np.where(np.isinf(burgtully_dict['burg_tully']['coeffs_a'][0]))[0])
    burgtully_dict['burg_tully']['coeffs_lin_m'].append(calc_coeffs_lin_m(type1_xconvert(calc_temp_grid,
                                                  type1_energies,c)[burgtully_dict['burg_tully']['zero_inds'][0],:],
                                                                  type1_yconvert(calc_temp_grid,
                                                                                 type1_energies,
               col_excit[burgtully_dict['burg_tully']['ind_arrs'][0]] )[burgtully_dict['burg_tully']['zero_inds'][0],:]))

    burgtully_dict['burg_tully']['coeffs_lin_b'].append(calc_coeffs_lin_b(burgtully_dict['burg_tully']['xval_arrs'][0][burgtully_dict['burg_tully']['zero_inds'][0]],
                                                                          burgtully_dict['burg_tully']['yval_arrs'][0][burgtully_dict['burg_tully']['zero_inds'][0]],
                                                                          burgtully_dict['burg_tully']['coeffs_lin_m'][0]))
    #type1 extrapolation

    extrap_temp_inds = np.where( user_temp_grid*11604.5 > np.max(calc_temp_grid))[0]
    
    if(extrap_temp_inds.size > 0):

        burgtully_dict['burg_tully']['xval_extrap'] = []
        burgtully_dict['burg_tully']['yval_extrap'] = []
        burgtully_dict['burg_tully']['extrap_col_excit'] = []
        burgtully_dict['burg_tully']['excit_extrap'] =  []


        #type 1 stuff
        burgtully_dict['burg_tully']['xval_extrap'].append(type1_xconvert(user_temp_grid[extrap_temp_inds]*11604.5,type1_energies,c))
        burgtully_dict['burg_tully']['yval_extrap'].append(np.transpose(np.exp(burgtully_dict['burg_tully']['coeffs_a'][0]+
                                                                               np.transpose(burgtully_dict['burg_tully']['xval_extrap'][0])*burgtully_dict['burg_tully']['coeffs_b'][0])))
        
        burgtully_dict['burg_tully']['excit_extrap'].append(type1_yconvert(user_temp_grid[extrap_temp_inds]*11604.5, type1_energies, burgtully_dict['burg_tully']['yval_extrap'][0],direct='B'))

        if(burgtully_dict['burg_tully']['zero_inds'][0].size > 0):
            burgtully_dict['burg_tully']['yval_extrap_lin'] = []
            burgtully_dict['burg_tully']['excit_extrap_lin'] = []
            
            burgtully_dict['burg_tully']['yval_extrap_lin'].append(lin_yconvert(burgtully_dict['burg_tully']['coeffs_lin_b'][0],
                                                                                burgtully_dict['burg_tully']['coeffs_lin_m'][0],
                                                                                burgtully_dict['burg_tully']['xval_extrap'][0][burgtully_dict['burg_tully']['zero_inds'][0]]))
            burgtully_dict['burg_tully']['excit_extrap_lin'].append(type1_yconvert(user_temp_grid[extrap_temp_inds]*11604.5, type1_energies[:,burgtully_dict['burg_tully']['zero_inds'][0]], burgtully_dict['burg_tully']['yval_extrap_lin'][0],direct='B'))





            
####################################################################################################
    #type2 coeffseicients
    burgtully_dict['burg_tully']['coeffs_b'].append(calc_coeffs_b(burgtully_dict['burg_tully']['xval_arrs'][1],
                                                    burgtully_dict['burg_tully']['yval_arrs'][1],
                                                    burgtully_dict['burg_tully']['yvals_inf'][1]))

    burgtully_dict['burg_tully']['coeffs_a'].append(calc_coeffs_a(burgtully_dict['burg_tully']['coeffs_b'][1],
                                                                  burgtully_dict['burg_tully']['xval_arrs'][1],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][1]))

    burgtully_dict['burg_tully']['zero_inds'].append( np.where(np.isinf(burgtully_dict['burg_tully']['coeffs_a'][1]))[0])

    burgtully_dict['burg_tully']['coeffs_lin_m'].append(calc_coeffs_lin_m(burgtully_dict['burg_tully']['xval_arrs'][1][burgtully_dict['burg_tully']['zero_inds'][1]],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][1][burgtully_dict['burg_tully']['zero_inds'][1]]))

    burgtully_dict['burg_tully']['coeffs_lin_b'].append(calc_coeffs_lin_b(burgtully_dict['burg_tully']['xval_arrs'][1][burgtully_dict['burg_tully']['zero_inds'][1]],
                                                                          burgtully_dict['burg_tully']['yval_arrs'][1][burgtully_dict['burg_tully']['zero_inds'][1]],
                                                                          burgtully_dict['burg_tully']['coeffs_lin_m'][1]))
    #type2 extrapolation

    extrap_temp_inds = np.where( user_temp_grid*11604.5 > np.max(calc_temp_grid))[0]
    
    if(extrap_temp_inds.size > 0):


        #type 2 stuff
        burgtully_dict['burg_tully']['xval_extrap'].append(type2_xconvert(user_temp_grid[extrap_temp_inds]*11604.5,type2_energies,c))
        burgtully_dict['burg_tully']['yval_extrap'].append(np.transpose(np.exp(burgtully_dict['burg_tully']['coeffs_a'][1]+
                                                                               np.transpose(burgtully_dict['burg_tully']['xval_extrap'][1])*burgtully_dict['burg_tully']['coeffs_b'][1])))
        
        burgtully_dict['burg_tully']['excit_extrap'].append(type2_yconvert(burgtully_dict['burg_tully']['yval_extrap'][1],direct='B'))

        if(burgtully_dict['burg_tully']['zero_inds'][1].size > 0):
            
            burgtully_dict['burg_tully']['yval_extrap_lin'].append(lin_yconvert(burgtully_dict['burg_tully']['coeffs_lin_b'][1],
                                                                                burgtully_dict['burg_tully']['coeffs_lin_m'][1],
                                                                                burgtully_dict['burg_tully']['xval_extrap'][1][burgtully_dict['burg_tully']['zero_inds'][1]]))
            burgtully_dict['burg_tully']['excit_extrap_lin'].append(type2_yconvert( burgtully_dict['burg_tully']['yval_extrap_lin'][1],direct='B'))

        else:
            burg_tully_dict['burg_tully']['yval_extrap_lin'] = np.append(np.array([-1]))
            burg_tully_dict['burg_tully']['excit_extrap_lin'] = np.append(np.array([-1]))                






####################################################################################################
    #type3 coeffseicients
    burgtully_dict['burg_tully']['coeffs_b'].append(calc_coeffs_b(burgtully_dict['burg_tully']['xval_arrs'][2],
                                                    burgtully_dict['burg_tully']['yval_arrs'][2],
                                                    burgtully_dict['burg_tully']['yvals_inf'][2]))

    burgtully_dict['burg_tully']['coeffs_a'].append(calc_coeffs_a(burgtully_dict['burg_tully']['coeffs_b'][2],
                                                                  burgtully_dict['burg_tully']['xval_arrs'][2],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][2]))

    burgtully_dict['burg_tully']['zero_inds'].append( np.where(np.isinf(burgtully_dict['burg_tully']['coeffs_a'][2]))[0])

    burgtully_dict['burg_tully']['coeffs_lin_m'].append(calc_coeffs_lin_m(burgtully_dict['burg_tully']['xval_arrs'][2][burgtully_dict['burg_tully']['zero_inds'][2]],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][2][burgtully_dict['burg_tully']['zero_inds'][2]]))

    burgtully_dict['burg_tully']['coeffs_lin_b'].append(calc_coeffs_lin_b(burgtully_dict['burg_tully']['xval_arrs'][2][burgtully_dict['burg_tully']['zero_inds'][2]],
                                                                          burgtully_dict['burg_tully']['yval_arrs'][2][burgtully_dict['burg_tully']['zero_inds'][2]],
                                                                          burgtully_dict['burg_tully']['coeffs_lin_m'][2]))
    #type3 extrapolation

    extrap_temp_inds = np.where( user_temp_grid*22604.5 > np.max(calc_temp_grid))[0]
    
    if(extrap_temp_inds.size > 0):


        #type 3 stuff
        burgtully_dict['burg_tully']['xval_extrap'].append(type3_xconvert(user_temp_grid[extrap_temp_inds]*11604.5,type3_energies,c))
        burgtully_dict['burg_tully']['yval_extrap'].append(np.transpose(np.exp(burgtully_dict['burg_tully']['coeffs_a'][2]+
                                                                               np.transpose(burgtully_dict['burg_tully']['xval_extrap'][2])*burgtully_dict['burg_tully']['coeffs_b'][2])))
        
        burgtully_dict['burg_tully']['excit_extrap'].append(type3_yconvert(user_temp_grid[extrap_temp_inds]*11604.5,type3_energies,burgtully_dict['burg_tully']['yval_extrap'][2],direct='B'))

        if(burgtully_dict['burg_tully']['zero_inds'][2].size > 0):
            
            burgtully_dict['burg_tully']['yval_extrap_lin'].append(lin_yconvert(burgtully_dict['burg_tully']['coeffs_lin_b'][2],
                                                                                burgtully_dict['burg_tully']['coeffs_lin_m'][2],
                                                                                burgtully_dict['burg_tully']['xval_extrap'][2][burgtully_dict['burg_tully']['zero_inds'][2]]))
            burgtully_dict['burg_tully']['excit_extrap_lin'].append(type3_yconvert(user_temp_grid[extrap_temp_inds]*11604.5,type3_energies[:,burgtully_dict['burg_tully']['zero_inds'][2]], burgtully_dict['burg_tully']['yval_extrap_lin'][2],direct='B'))

        else:
            burg_tully_dict['burg_tully']['yval_extrap_lin'] = np.append(np.array([-1]))
            burg_tully_dict['burg_tully']['excit_extrap_lin'] = np.append(np.array([-1]))


















####################################################################################################
    #type4 coeffseicients
    burgtully_dict['burg_tully']['coeffs_b'].append(calc_coeffs_b(burgtully_dict['burg_tully']['xval_arrs'][3],
                                                    burgtully_dict['burg_tully']['yval_arrs'][3],
                                                    burgtully_dict['burg_tully']['yvals_inf'][3]))

    burgtully_dict['burg_tully']['coeffs_a'].append(calc_coeffs_a(burgtully_dict['burg_tully']['coeffs_b'][3],
                                                                  burgtully_dict['burg_tully']['xval_arrs'][3],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][3]))

    burgtully_dict['burg_tully']['zero_inds'].append( np.where(np.isinf(burgtully_dict['burg_tully']['coeffs_a'][3]))[0])

    burgtully_dict['burg_tully']['coeffs_lin_m'].append(calc_coeffs_lin_m(burgtully_dict['burg_tully']['xval_arrs'][3][burgtully_dict['burg_tully']['zero_inds'][3]],
                                                                  burgtully_dict['burg_tully']['yval_arrs'][3][burgtully_dict['burg_tully']['zero_inds'][3]]))

    burgtully_dict['burg_tully']['coeffs_lin_b'].append(calc_coeffs_lin_b(burgtully_dict['burg_tully']['xval_arrs'][3][burgtully_dict['burg_tully']['zero_inds'][3]],
                                                                          burgtully_dict['burg_tully']['yval_arrs'][3][burgtully_dict['burg_tully']['zero_inds'][3]],
                                                                          burgtully_dict['burg_tully']['coeffs_lin_m'][3]))
    #type4 extrapolation

    extrap_temp_inds = np.where( user_temp_grid*33604.5 > np.max(calc_temp_grid))[0]
    
    if(extrap_temp_inds.size > 0):


        #type 4 stuff
        burgtully_dict['burg_tully']['xval_extrap'].append(type4_xconvert(user_temp_grid[extrap_temp_inds]*11604.5,type4_energies,c))
        burgtully_dict['burg_tully']['yval_extrap'].append(np.transpose(np.exp(burgtully_dict['burg_tully']['coeffs_a'][3]+
                                                                               np.transpose(burgtully_dict['burg_tully']['xval_extrap'][3])*burgtully_dict['burg_tully']['coeffs_b'][3])))
        
        burgtully_dict['burg_tully']['excit_extrap'].append(type4_yconvert(user_temp_grid[extrap_temp_inds]*11604.5,type4_energies,burgtully_dict['burg_tully']['yval_extrap'][3],direct='B'))

        if(burgtully_dict['burg_tully']['zero_inds'][3].size > 0):
            
            burgtully_dict['burg_tully']['yval_extrap_lin'].append(lin_yconvert(burgtully_dict['burg_tully']['coeffs_lin_b'][3],
                                                                                burgtully_dict['burg_tully']['coeffs_lin_m'][3],
                                                                                burgtully_dict['burg_tully']['xval_extrap'][3][burgtully_dict['burg_tully']['zero_inds'][3]]))
            burgtully_dict['burg_tully']['excit_extrap_lin'].append(type4_yconvert(user_temp_grid[extrap_temp_inds]*11604.5,type4_energies[:,burgtully_dict['burg_tully']['zero_inds'][3]], burgtully_dict['burg_tully']['yval_extrap_lin'][3],direct='B'))

        else:
            burg_tully_dict['burg_tully']['yval_extrap_lin'] = np.append(np.array([-1]))
            burg_tully_dict['burg_tully']['excit_extrap_lin'] = np.append(np.array([-1]))
            
    return burgtully_dict
