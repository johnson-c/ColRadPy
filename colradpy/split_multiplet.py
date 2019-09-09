import numpy  as np

def split_condon_shortley(s1,l1,j1,s2,l2,j2):
    """ This function calculates splitting percentages for LS
        coupled PECs that are to be split to LSJ.
        This is based on the orginal 'split_multiplet.pro' from ADAS, orginal
        author Martin O'Mulane. Page 238 The Theory of Atomic Spectra, Condon, Shortley 1935
        
    Args:
      :param s1: 'S' quantum number array for the PECs from the lower term
      :type s1: integer

      :param l1: 'L' quantum number array for the PECs from the lower term
      :type l1: integer

      :param s2: 'S' quantum number array for the PECs from the upper term
      :type s2: integer

      :param l2: 'L' quantum number array for the PECs from the upper term
      :type l2: integer

    :returns: integer splitting percentage for the term resolved PEC
    """
    
    if( (j1 == j2-1) and (l1 == l2+1)):
        res = (j1+s1-l1+1) * (l1+s1-j1) * (j1+s1-l1+2)  * (l1+s1-j1-1) / (4.0 * (j1+1) )
    elif((j1 == j2)   and (l1 == l2+1)):
        res = (2*j1+1) * (j1+l1-s1) * (j1+s1-l1+1) * (s1+l1+1+j1) * (s1+l1-j1) / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2+1)):
        res = (j1+l1-s1-1) * (j1+l1-s1) * (s1+l1+j1+1) * (s1+l1+j1) / ( 4.0*j1 )

    elif((j1 == j2-1) and (l1 == l2)):
        res = (j1-s1+l1+1) * (j1+s1-l1+1) * (s1+l1+j1+2) * (s1+l1-j1) / ( 4.0*(j1+1) )
    elif((j1 == j2)   and (l1 == l2)):
        res = (2*j1+1) * ( j1*(j1+1) - s1*(s1+1) + l1*(l1+1) )**2.0 / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2)):
        res = (j1-s1+l1) * (j1+s1-l1) * (s1+l1+j1+1) * (s1+l1+1-j1) / ( 4.0*j1 )
        
    elif((j1 == j2-1) and (l1 == l2-1)):
        res = (j1+l1-s1-1) * (j1+l1-s1) * (s1+l1+j1+1) * (s1+l1+j1) / ( 4.0*j1 )
    elif((j1 == j2)   and (l1 == l2-1)):
        res = (2*j1+1) * (j1+l1-s1) * (j1+s1-l1+1) * (s1+l1+1+j1) * (s1+l1-j1) / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2-1)):
        res = (j1+s1-l1+1) * (l1+s1-j1) * (j1+s1-l1+2)  * (l1+s1-j1-1) / (4.0 * (j1+1) )
    else:
        res = []

    return res



def split_multiplet(s_low,l_low,s_up,l_up):

    """ This function calculates splitting percentages for LS
        coupled PECs that are to be split to LSJ.
        This is based on the orginal 'split_multiplet.pro' from ADAS, orginal
        author Martin O'Mulane. 

    Args:
      :param s_low: 'S' quantum number array for the PECs from the lower term
      :type s_low: integer array

      :param l_low: 'L' quantum number array for the PECs from the lower term
      :type l_low: integer array

      :param s_up: 'S' quantum number array for the PECs from the upper term
      :type s_up: integer array

      :param l_up: 'L' quantum number array for the PECs from the upper term
      :type l_up: integer array

    :returns: array 'j' value array for the upper level
              array 'j' value array for the lower level
              array splitting percentages for the term resolved PEC
    """
    j_up = []
    j_low = []
    res = []
    #finding all possible j values for upper and lower

    #l_low has to be smaller than l_up
    switch = False
    if(l_up < l_low):
        l_t = l_low
        s_t = s_low
        s_low = s_up
        l_low = l_up
        l_up = l_t
        s_up = s_t
        switch=True
    j1_arr = np.arange(np.abs(l_up-s_up),np.abs(l_up+s_up)+1.,1.)
    j2_arr = np.arange(np.abs(l_low-s_low),np.abs(l_low+s_low)+1.,1.)

    for j1 in j1_arr:
        for j2 in j2_arr:
            j_up.append(j1)
            j_low.append(j2)
            if( (np.abs(j2-j1) < 2) and (np.abs(l_up - l_low) < 2) and (np.abs(j1+j2) >= 1)):
                res.append( split_condon_shortley(s_up,l_up,j1,s_low,l_low,j2))
            else:
                res.append(0)


    if(switch):
        return np.asarray(j_low),np.asarray(j_up),np.asarray(res)        
    else:
        return np.asarray(j_up),np.asarray(j_low),np.asarray(res)
