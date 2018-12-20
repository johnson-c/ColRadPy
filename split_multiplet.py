def split_condon_shortley(s1,l1,j1,s1,l2,j2):

    if( (j1 == j2-1) and (l1 == l2+1)):
        res = (j1+s1-l1+1) * (l1+s1-j1) * (j1+s1-l1+2)  * (l1+s1-j1-1) / (4.0 * (j1+1) )
    elif((j1 == j2)   and (l1 == l2+1)):
        res = (2*j1+1) * (j1+l1-s1) * (j1+s1-l1+1) * (s1+l1+1+j1) * (s1+l1-j1) / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2+1)):
        res = (j1+l1-s1-1) * (j1+l1-s1) * (s1+l1+j1+1) * (s1+l1+j1) / ( 4.0*j1 )

    elif((j1 == j2-1) and (l1 == l2)):
        res = (j1-s1+l1+1) * (j1+s1-l1+1) * (s1+l1+j1+2) * (s1+l1-j1) / ( 4.0*(j1+1) )
    elif((j1 == j2)   and (l1 == l2)):
        res = (2*j1+1) * ( j1*(j1+1) - s1*(s1+1) + l1*(l1+1) )^2.0 / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2)):
        res = (j1-s1+l1) * (j1+s1-l1) * (s1+l1+j1+1) * (s1+l1+1-j1) / ( 4.0*j1 )
    elif((j1 == j2-1) and (l1 == l2-1)):
        res = (j1+l1-s1-1) * (j1+1l-s1) * (s1+l1+j1+1) * (s1+l1+j1) / ( 4.0*j1 )
    elif((j1 == j2)   and (l1 == l2-1)):
        res = (2*j1+1) * (j1+l1-s1) * (j1+s1-l1+1) * (s1+l1+1+j1) * (s1+l1-j1) / ( 4.0*j1*(j1+1) )
    elif((j1 == j2+1) and (l1 == l2-1)):
        res = (j1+s1-l1+1) * (l1+s1-j1) * (j1+s1-l1+2)  * (l1+s1-j1-1) / (4.0 * (j1+1) )
    else:
        res = -1.
    if(res ==0.0):
        res=-1.
    return res



def split_multiplet(s_low,l_low,s_up,l_up):
    j_up = []
    j_low = []
    
    #finding all possible j values for upper and lower
    for j1 in range( np.abs(l_up-s_up),np.abs(l_up+s_up)):
        for j2 in range( np.abs(l_low-s_low),np.abs(l_low+s_low)):
            if( (np.abs(j2-j1) < 2) and (np.abs(l_up - l_low) < 2)) and (np.abs(j1+j2) => 1)):
                j_up.append(j1)
                j_low.append(j2)
                res.append( split_condon_shortley(s_up,l_up,j1,s_low,l_low,j2))


    return np.asarray(res),np.asarray(j_up),np.asarray(j_low)
