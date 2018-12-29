import numpy  as np

def split_condon_shortley(s1,l1,j1,s2,l2,j2):

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
    if(res ==0.0):
        res=[]
    return res



def split_multiplet(s_low,l_low,s_up,l_up):
    j_up = []
    j_low = []
    res = []
    #finding all possible j values for upper and lower

    #l_low has to be smaller than l_up
    if(l_up < l_low):
        l_t = l_low
        s_t = s_low
        s_low = s_up
        l_low = l_up
        l_up = l_t
        s_up = s_t

    j1_arr = np.arange(np.abs(l_up-s_up),np.abs(l_up+s_up)+1.,1.)
    j2_arr = np.arange(np.abs(l_low-s_low),np.abs(l_low+s_low)+1.,1.)
    
    for j1 in j1_arr:
        for j2 in j2_arr:

            if( (np.abs(j2-j1) < 2) and (np.abs(l_up - l_low) < 2) and (np.abs(j1+j2) >= 1)):
                j_up.append(j1)
                j_low.append(j2)
                res.append( split_condon_shortley(s_up,l_up,j1,s_low,l_low,j2))
            

    return np.asarray(j_up),np.asarray(j_low),np.asarray(res)



'''

he2 = colradpy('/home/curtis/adf04_files/mom97_ls#he1.dat',[0],np.array([100]),np.array([1.e13]),use_recombination=False,use_recombination_three_body = False,use_ionization=True)
he2.solve_cr()
config_l = []
config_u = []    
L_l = []
L_u = []
S_l = []
S_u = []
j_l = []
j_u = []
ress = []
pecs = []
for i in range(0,len(he2.data['processed']['pec_levels'])):
    up = he2.data['processed']['pec_levels'][i,0]
    low = he2.data['processed']['pec_levels'][i,1]
    ju,jl,res = split_multiplet( (he2.data['atomic']['S'][up]-1)/2.,he2.data['atomic']['L'][up],(he2.data['atomic']['S'][low]-1)/2.,he2.data['atomic']['L'][low])
    j_l.append(jl)
    j_u.append(ju)
    ress.append(res)
    if(res.size>0):
        if( -1 not in res):
            pecs.append(np.einsum('ijk,l->lijk',he2.data['processed']['pecs'][i],res/np.sum(res)))
        else:
            pecs.append(he2.data['processed']['pecs'][i])
    else:
        pecs.append(he2.data['processed']['pecs'][i])    

    ress.append(res)
    if(res.size>0):
        if(-1 not in res):
            for j in range(0,len(ju)):
                j_l.append(jl[j])
                j_u.append(ju[j])
                pecs.append(he2.data['processed']['pecs'][i]*res[j]/np.sum(res))
                config_u.append(he2.data['atomic']['config'][up])
                config_l.append(he2.data['atomic']['config'][low])
                L_l.append(he2.data['atomic']['L'][low])
                L_u.append(he2.data['atomic']['L'][up])
                S_u.append(he2.data['atomic']['S'][up])
                S_l.append(he2.data['atomic']['S'][low])
        else:
                j_l.append(jl[0])
                j_u.append(ju[0])
                pecs.append(he2.data['processed']['pecs'][i])
                config_u.append(he2.data['atomic']['config'][up])
                config_l.append(he2.data['atomic']['config'][low])
                L_l.append(he2.data['atomic']['L'][low])
                L_u.append(he2.data['atomic']['L'][up])
                S_u.append(he2.data['atomic']['S'][up])
                S_l.append(he2.data['atomic']['S'][low])
    else:
        j_l.append(he2.data['atomic']['w'][up])
        j_u.append(he2.data['atomic']['w'][low])
        pecs.append(he2.data['processed']['pecs'][i])        
        config_u.append(he2.data['atomic']['config'][up])
        config_l.append(he2.data['atomic']['config'][low])
        L_l.append(he2.data['atomic']['L'][low])
        L_u.append(he2.data['atomic']['L'][up])
        S_u.append(he2.data['atomic']['S'][up])
        S_l.append(he2.data['atomic']['S'][low])
        
pecs = np.asarray(pecs)
    '''
