################################################################################
# file name         : read_adf04.py
# author            : Curt Johnson
# description       : This code reads adf04 formatted files
# version           : 0.1
# python version    : 2.7.12 ipython 2.4.1
# dependencies      : numpy, re
# license           : to be freed after code is stable and papers are published
#
# todo          check with all the adf04 files in adas
#
#
################################################################################

import numpy as np
import re
from operator import xor

def read_adf04(file_path):
    file_name = file_path
    dict = {}
    f = open(file_name)

    ############################################################################
    #
    # reading the first line
    #
    ############################################################################
    first_line = re.split('[ ]', f.readline())

    ion_pot = []
    ion_term = []
    for i in range(0,len(first_line)):

        if(first_line[i] != '' and  not dict.has_key('element')):
            dict['element'] = first_line[i]
            if(len( dict['element']) > 3):
                tt = re.split('[+]',dict['element'])[0]
                dict['element'] = tt[0]
                dict['charge_state'] = tt[1]
        elif(first_line[i] != '' and  dict.has_key('element') and \
             not dict.has_key('charge_state') and first_line[i] !='+'):
            dict['charge_state'] = first_line[i]

        elif(first_line[i] != '' and  dict.has_key('element') and \
           dict.has_key('charge_state') and not dict.has_key('iz0')):

            dict['iz0'] = first_line[i]
        elif(first_line[i] != '' and  dict.has_key('element') and \
           dict.has_key('charge_state') and dict.has_key('iz0') and \
           not dict.has_key('iz1')):
            dict['iz1'] = first_line[i]
        elif(first_line[i] != '' and  dict.has_key('element') and \
           dict.has_key('charge_state') and dict.has_key('iz0') and \
           dict.has_key('iz1') and not dict.has_key('ion_term')):
             ion_pot.append(float(re.split('[(]',first_line[i])[0]))
             ion_term.append(re.split('[)]',re.split('[(]',first_line[i])[1])[0])
    dict['ion_pot'] = np.asarray(ion_pot)
    dict['ion_term'] = np.asarray(ion_term)
    ############################################################################
    #
    # done reading the first line, reading the level information
    #
    ############################################################################

    dict['charge_state'] = int(dict['charge_state'])
    config = []
    L = []
    S = []
    w = []
    energy = []
    zpla = []
    zpla1 = []
    zplaa = []
    zp_arr = np.ones(len(dict['ion_term']))
    zp_arr = zp_arr*-1.
    neg_found = False
    level_ind_found = False
    conf_temp = ''
    conf_tmp_found = False
    maybe_w=''
    ttttmp = f.readline()
    level_line = re.split('[ ]',ttttmp)

    while(  xor( ('-1\n' not in level_line),(level_line[0] == 'C') )):
        neg_found = False
        level_ind_found = False
        conf_temp = ''
        conf_tmp_found = False
        maybe_w=''
        energy_found = False
        zplaf = False
        conf_tmp = ''

        for i in range(0,len(level_line)):

            if(level_line[i] != '' and  not level_ind_found ):
                level_ind_found = True
            elif(level_line[i] != '' and  level_ind_found and\
                 level_line[i][0] != '(' and not conf_tmp_found):

                conf_tmp = conf_tmp + level_line[i] + ' '
            elif(level_line[i] != '' and  level_ind_found and\
                 level_line[i][0] == '('):
                conf_tmp_found = True
                sl = re.split('[(]',level_line[i])[1]
                S.append( int(re.split('[)]',sl)[1]))
                L.append( int(re.split('[)]',sl)[0]))
                maybe_w = re.split('[(]',level_line[i])[-1]
                if(maybe_w):
                    w.append(float(re.split('[)]',maybe_w)[0]))
            elif(not maybe_w and conf_tmp_found and level_ind_found
                 and level_line[i] != ''):
                if(level_line[i][-1] == ')'):
                    w.append(float(re.split('[)]',level_line[i])[0]))
                    maybe_w = 'found'
            elif( maybe_w and conf_tmp_found and level_ind_found and
                  level_line[i]!='' and level_line[i][0]!='{' and not energy_found):
                energy.append(float(level_line[i]))
                energy_found = True
                zpla1.append(-1)
            elif( maybe_w and conf_tmp_found and level_ind_found and
                  level_line[i]!='' and level_line[i][0] == '{'):
                if(level_line[i][1] !='X'):
                    ion_sig = int(level_line[i][1])-1
                    zp_arr[ion_sig] = float(re.split('[}]',level_line[i])[1])

                if(level_line[i][1]=='1'):
                    zpla.append( float(re.split('[}]',level_line[i])[1]))
                    zplaf=True
                    
                else:
                    zpla.append(-1)
                    zplaf=True
            
            if(maybe_w and conf_tmp_found and level_ind_found and
                  level_line[i]!='' and level_line[i][0] == '{' and level_line[i][1]=='2'):
                    zpla1[len(energy)-1] =  float(re.split('[}]',level_line[i])[1])

                    
        
        if( '\n' in level_line[i]):

            zplaa.append(zp_arr)
            zp_arr = np.ones(len(dict['ion_term']))
            zp_arr = zp_arr*-1.        
        config.append(conf_tmp[0:len(conf_tmp)-1])
        ttttmp = f.readline()        
        level_line = re.split('[ ]',ttttmp)
    dict['zpla'] = np.asarray(zplaa)
    ############################################################################
    #
    # end reading level info start reading temperature grid
    #
    ############################################################################
    t = f.readline()
    temp_grid_tmp = re.split('[ ]',t)
    one = False
    two = False
    adf04_temp_grid = []
    for i in range(0,len(temp_grid_tmp)):

        if(temp_grid_tmp[i] !='' and one and two and temp_grid_tmp[i] !='\n'):
            if( '+' in temp_grid_tmp[i]):
                temp_grid_tmp[i] = temp_grid_tmp[i].replace('+','E')
            elif('-' in temp_grid_tmp[i]):
                temp_grid_tmp[i] = temp_grid_tmp[i].replace('-','E-')
            adf04_temp_grid.append(float(temp_grid_tmp[i]))

        if(temp_grid_tmp[i] !='' and one and not two):
            two = True

        if(temp_grid_tmp[i] !='' and not one):
            one = True

    dict['temp_grid'] = np.asarray(adf04_temp_grid)
    dict['energy'] =np.asarray(energy)
    dict['L'] = np.asarray(S)
    dict['S'] = np.asarray(L)
    dict['w'] = np.asarray(w)
    dict['config'] = np.asarray(config)

    dict['zpla1'] = np.asarray(zpla1)
    ################################################################################
    #
    # end reading temp grid start reading in the atomic data from the adf04 file
    #
    ################################################################################
    a = np.loadtxt(file_name,skiprows=len(L)+3,dtype=str,delimiter='\n')

    #declare the colisional stuff
    col_excit = []
    col_transitions = []
    col_trans = np.array([-1,-1])
    #declare the recombination stuff
    recomb_excit = []
    recomb_transitions = []
    recomb_trans = np.array([-1,-1])
    #declare the ion stuff
    ion_excit = []
    ion_transitions = []
    ion_trans = np.array([-1,-1])
    #declare the a value stuff
    a_val = []
    inf_engy= []

    i = 0 #this is the index of the line in the adf04 file starting at the colision
    b = re.split('[ ]',a[i])#this will the the line that gets read from the loadtxt
    while( '-1' not in b):

        #################################################
        #
        # get the excitation and a values from the adf04
        #
        ##################################################

        first_ind = -1  #tells if first index found
        second_ind = -1 #tells if second index found
        if(b[0] == ''): #this means that it is collisional excitation and A vals
            k = 0
            #going to loop over to keep looping until the indices are found
            #this is happending one ' ' at a time
            while(second_ind != 1): #for k in range(1,len(b)):
                #check for the second index first 
                if(b[k] !='' and first_ind==1 and second_ind == 0):

                    col_trans[1] = int(b[k])
                    second_ind=1
                #check if the first index found
                if(b[k] !='' and first_ind ==-1):
                    col_trans[0] = int(b[k])
                    first_ind = 1
                    second_ind=0
                k = k + 1 # that ' ' step
            col_transitions.append(np.array([col_trans[0],col_trans[1]]))#append the indexes
            problem_w_b = np.where(np.asarray(b[k:]) =='')[0]
            if(len(problem_w_b) > 0):
                b = b[k:problem_w_b[0]]
            else:
                b = b[k:] #now b is just the aval and the upsilons
            c = np.zeros(len(b)) #make a matrix for the upsilons
            if(len(b) <len(adf04_temp_grid)+2):

                tmpp = b[-1][7:]
                b[-1] = b[-1][0:7]
                b.append(tmpp)


            #this is because adas has poor scientific notation
            for j in range(0,len(adf04_temp_grid)+1): # +1 because the a vale in there too
                if( '+' in b[j]):
                    b[j] = b[j].replace('+','E')
                elif('-' in b[j]):
                    b[j] = b[j].replace('-','E-')

                c[j] =float(b[j])
            a_val.append(c[0])#append the A values
            col_excit.append(c[1:1+len(adf04_temp_grid)])#append the upsilon matrix


            b = np.asarray(b)
            if(np.where(b == '')[0]):
                b = np.delete(b,np.where(b == '')[0])

            if(len(b) == len(adf04_temp_grid)+2):
                inf_e = b[-1]
                l = []
                l.append(inf_e[0])
                if( '+' in inf_e[1:]):
                    l.append(inf_e[1:].replace('+','E'))
                elif('-' in inf_e[1:]):
                    l.append(inf_e[1:].replace('-','E-'))

                inf_engy.append(float(''.join(l)))
                c[j] =float(b[j])
        ######################################################################
        #
        # end excitation and a values from the adf04, start looking for recomb
        #
        ######################################################################

        elif(b[0] == 'R'):
            k = 0
            recomb_first_ind =  -1
            recomb_second_ind = -1
            while(recomb_second_ind != 1):
                if(b[k] != 'R' and b[k] != '' and recomb_first_ind==1 and recomb_second_ind==0):

                    recomb_trans[0] = int(b[k])
                    recomb_second_ind = 1
                if(b[k] != 'R' and b[k] != '' and recomb_first_ind==-1):
                    recomb_trans[1] = int(b[k])
                    recomb_first_ind = 1
                    recomb_second_ind = 0
                k = k + 1 #that ' ' step
            recomb_transitions.append(np.array([recomb_trans[0],recomb_trans[1]]))

            b = b[k:]
            #get rid of '' leading up to the upsilons
            upsilon_found = False
            ups_ind = 0 
            while( not upsilon_found):
                if(b[ups_ind] != ''):
                    upsilon_found = True
                ups_ind = ups_ind + 1
            b = b[ups_ind-1:]

            c = np.zeros(len(b)) #make a matrix for the upsilons
            #this is because adas has poor scientific notation
            for j in range(0,len(b)):
                if( '+' in b[j]):
                    b[j] = b[j].replace('+','E')
                elif('-' in b[j]):
                    b[j] = b[j].replace('-','E-')
                c[j] =float(b[j])
            recomb_excit.append(c)


        elif(b[0] == 'S'):
            k = 0
            ion_first_ind =  -1
            ion_second_ind = -1
            while(ion_second_ind != 1):
                if(b[k] != 'R' and b[k] != '' and ion_first_ind==1 and ion_second_ind==0):

                    ion_trans[0] = int(b[k])
                    ion_second_ind = 1
                if(b[k] != 'S' and b[k] != '' and ion_first_ind==-1):
                    ion_trans[1] = int(b[k])
                    ion_first_ind = 1
                    ion_second_ind = 0
                k = k + 1 #that ' ' step

            ion_transitions.append(np.array([ion_trans[1],ion_trans[0]]))

            b = b[k:]
            #get rid of '' leading up to the upsilons
            upsilon_found = False
            ups_ind = 0 
            while( not upsilon_found):
                if(b[ups_ind] != ''):
                    upsilon_found = True
                ups_ind = ups_ind + 1
            b = b[ups_ind-1:]

            c = np.zeros(len(b)) #make a matrix for the upsilons
            #this is because adas has poor scientific notation
            for j in range(0,len(b)):
                if( '+' in b[j]):
                    b[j] = b[j].replace('+','E')
                    c[j] =float(b[j])
                elif('-' in b[j]):
                    b[j] = b[j].replace('-','E-')
                    c[j] =float(b[j])
                    
            ion_excit.append(c)
        i = i + 1 #incrimententing the next line in the file index
        b = re.split('[ ]',a[i])#getting the next line in the file done this way
        #so that we can find '-1' signifying the end of the file before we loop




    a_val = np.asarray(a_val)#make the a vals a np array instead of list
    #make the upsilons np array instead of list
    dict['a_val'] = a_val
    dict['col_transitions'] = np.asarray(col_transitions)
    dict['recomb_transitions'] = np.asarray(recomb_transitions)
    dict['ion_transitions'] = np.asarray(ion_transitions)
    dict['col_excit'] = np.asarray(col_excit)
    dict['recomb_excit'] = np.asarray(recomb_excit)
    dict['ion_excit'] = np.asarray(ion_excit)
    dict['inf_engy'] = np.asarray(inf_engy)
    return dict
