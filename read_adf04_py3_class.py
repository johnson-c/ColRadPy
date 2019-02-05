################################################################################
# file name         : read_adf04_py3.py
# author            : Curt Johnson
# description       : This code reads adf04 formatted files
# version           : 0.9
# python version    : 3.6.5 |Anaconda, INC.| ipython 6.4.0
# dependencies      : numpy, re, sys
# license           : to be freed after code is stable and papers are published
#
# todo          check with all the adf04 files in adas
#
#
################################################################################

import numpy as np
import re
import sys
import copy


def first_substring(strings, substring):
    return next(i for i, string in enumerate(strings) if substring in string[0])

def config_format(config_arr):
    conf = config_arr[0]
    if(len(config_arr)>1):
        for i in range(1,len(config_arr)):
            conf = conf+'.' + config_arr[i]
    return conf

def read_adf04(fil):
    file_path = fil

    file_name = file_path
    adf04 = {}
    adf04['atomic'] = {}
    adf04['input_file'] = {}
    adf04['rates'] = {}
    f = open(file_name)
    ############################################################################
    #
    # reading the first line
    #
    ############################################################################
    try:
        first_line = f.readline()
        adf04['atomic']['element'] = re.split('[+]',first_line)[0]
        first_linep = np.asarray(re.split('[ ]',re.split('[+]',first_line)[1]))
        first_linep_inds = np.where(first_linep !='')[0]
        adf04['atomic']['charge_state'] = int(first_linep[first_linep_inds[0]])
        adf04['atomic']['iz0'] = int(first_linep[first_linep_inds[1]])
        adf04['atomic']['iz1'] = int(first_linep[first_linep_inds[2]])
        adf04['atomic']['ion_pot'] = np.zeros(len(first_linep_inds)-3)

        ion_term = []
        for i in range(0,len(adf04['atomic']['ion_pot'])):
            adf04['atomic']['ion_pot'][i] =float(
                re.split('[(]',first_linep[first_linep_inds[i+3]])[0])
            ion_term.append(re.split(
                '[)]',re.split('[(]',first_linep[first_linep_inds[i+3]])[1])[0])

        adf04['atomic']['ion_term'] = np.asarray(ion_term)

    except:
        print("**Failed first line of the adf04**\n\n"+first_line+"\n\nfile is formatted incorrectly failed line.\nCHECK WHITE SPACE.\nIt should be formatted as below (Neutral Tungsten example). \nIf more than one ionization metastable,\nadd more ionization potentials and terms after the first\nW+ 0        74         1       63427.7(6D3)")
        sys.exit(1)

    ############################################################################
    #
    # done reading the first line, reading the level information
    #
    ############################################################################
    try:
        config = []
        L = []
        S = []
        w = []
        energy = []
        zpla = []
        zpla1 = []
        level_num = 0
        tmp_line =   first_line
        zpla = []
        zpla1 = []
        zz_check = np.array([-1])
        while('-1' not in tmp_line or '-1\n' not in tmp_line):
            level_num = level_num + 1
            tmp_line = f.readline()
            tmp = np.asarray(re.split('[  ]',tmp_line))

            if('-1' in tmp or '-1\n' in tmp):
                break
            else:
                tmp_inds = np.where(tmp !='')[0]
                #finding end of config
                config_stop = first_substring(tmp[tmp_inds],'(')

                config.append(config_format(tmp[tmp_inds[1:config_stop]]))
                S.append(int(re.split('[)]',re.split('[(]',tmp[tmp_inds[config_stop]])[1])[0]))
                L.append(int(re.split('[)]',re.split('[(]',tmp[tmp_inds[config_stop]])[1])[1]))

                if( tmp[tmp_inds[config_stop]].count(')') > 1):
                    w.append(float(re.split('[(]',tmp[tmp_inds[config_stop]])[2][0:4]))
                    offset = -1
                else:
                    w.append(float(re.split('[)]',re.split('[(]',tmp[tmp_inds[config_stop+1]])[0])[0]))
                    offset=0
                    
                energy.append(float(tmp[tmp_inds[config_stop+2+offset]]))
                zz = np.ones(len(adf04['atomic']['ion_pot']))*-1
                zz1 = np.ones(len(adf04['atomic']['ion_pot']))*-1
                

                for i in range(config_stop+3+offset,len(tmp_inds)):
                    if ('X' not in tmp[tmp_inds[i]]):
                        zpla1_line = int(re.split('[}]',tmp[tmp_inds[i]])[0][1])
                        zpla_line = float(re.split('[}]',tmp[tmp_inds[i]])[1])
                        if(zpla1_line not in zz_check):
                            if(zz_check[0] < 0):
                                zz_check[0] = int(zpla1_line)
                            else:
                                zz_check = np.append(zz_check,zpla1_line)

                        
                        zz[np.where(zz_check==zpla1_line)[0]] = zpla1_line
                        zz1[np.where(zz_check==zpla1_line)[0]] = zpla_line

                zpla.append(zz1)
                zpla1.append(zz)

        adf04['atomic']['config'] = np.asarray(config)
        adf04['atomic']['L'] = np.asarray(L)
        adf04['atomic']['S'] = np.asarray(S)
        adf04['atomic']['w'] = np.asarray(w)
        adf04['atomic']['energy'] = np.asarray(energy)
        adf04['atomic']['zpla'] = np.asarray(zpla)
        adf04['atomic']['zpla1'] = np.asarray(zpla1)
        adf04['atomic']['ion_pot_lvl'] = zz_check
    except:
        print("\n**FAILED LEVEL " +str(level_num)+"**\n\n"+"\n\nOne or more configuration/energy lines in the adf04 file are formatted incorrectly.\nIt should be formatted as below (C+3 example). \nIf more than one ionization metastable,\nadd more apla,zpla after the first ie. {2}1.00\n      2 1S2 2P1           (2)1( 2.5)    64555.4   {1}1.000")
        sys.exit(1)
    ############################################################################
    #
    # end reading level info start reading temperature grid
    #
    ############################################################################
    try:
        temp_line = f.readline()
        temp_line_arr = np.asarray(re.split('[ ]',temp_line))
        temp_line_inds = np.where( temp_line_arr!='')[0]
        adf04_temp_grid = []
        for i in range(2,len(temp_line_inds)):
            if( '+' in temp_line_arr[temp_line_inds[i]]):
                temp_tmp = temp_line_arr[temp_line_inds[i]].replace('+','E')
            elif('-' in temp_line_arr[temp_grid_inds[i]]):
                temp_tmp = temp_line_arr[temp_line_inds[i]].replace('-','E-')
            adf04_temp_grid.append(float(temp_tmp))
        adf04['input_file']['temp_grid'] = {}
        adf04['input_file']['temp_grid'] = np.asarray(adf04_temp_grid)
    except:
        print("\n**FAILED TEMPERATURE GRID**\n\n" +temp_line+"\n\nadf04 temperature grid formatted incorrectly.\nShould be formatted as below (C+3 example).\n  4.0    3       8.00+03 1.60+04 3.20+04 8.00+04 1.60+05 3.20+05 8.00+05 1.60+06 3.20+06 8.00+06 1.60+07")
        sys.exit(1)
    ################################################################################
    #
    # end reading temp grid start reading in the atomic data from the adf04 file
    #
    ################################################################################
    #declare the colisional stuff
    col_excit = []
    col_transitions = []
    #declare the recombination stuff
    recomb_excit = []
    recomb_transitions = []
    #declare the ion stuff
    ion_excit = []
    ion_transitions = []
    #declare the a value stuff
    a_val = []
    inf_engy= []
    while('-1' not in tmp or '-1\n' not in tmp):
        level_num = level_num + 1
        tmp_line = f.readline()
        tmp = np.asarray(re.split('[  ]',tmp_line))
        if('-1' in tmp or '-1\n' in tmp):
            break
        if(tmp[0] ==''):
            tmp_inds = np.where(tmp!='')[0]
            # Transition numbers
            col_transitions.append(np.array([int(tmp[tmp_inds[0]]),int(tmp[tmp_inds[1]])]))
            # A-values
            if( '+' in tmp[tmp_inds[2]]):
                a_val.append(float(tmp[tmp_inds[2]].replace('+','E')))
            elif('-' in tmp[tmp_inds[2]]):
                a_val.append(float(tmp[tmp_inds[2]].replace('-','E-')))
            #collisional values
            col_excit_row = np.zeros(len(adf04['input_file']['temp_grid']))
            for i in range(0,len(tmp_inds) - 3):

                if(len(tmp[tmp_inds[i+3]]) < 9):

                    if(i>len(adf04['input_file']['temp_grid'])-1):
                        if( '+' in tmp[tmp_inds[i+3]]):
                            inf_engy.append(float(tmp[tmp_inds[i+3]].replace('+','E')))
                        elif('-' in tmp[tmp_inds[i+3]]):
                            inf_engy.append(float(tmp[tmp_inds[i+3]].replace('-','E-')))
                    else:
                        if( '+' in tmp[tmp_inds[i+3]]):
                            col_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('+','E'))
                        elif('-' in tmp[tmp_inds[i+3]]):
                            col_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('-','E-'))

                        
                else:

                    last_temp = tmp[tmp_inds[i+3]][0:7]
                    inf_temp =  tmp[tmp_inds[i+3]][8:15]
                    if( '+' in last_temp):
                        col_excit_row[i] = float(last_temp.replace('+','E'))
                    elif('-' in last_temp):
                        col_excit_row[i] = float(last_temp.replace('-','E-'))
                    
                    if( '+' in inf_temp):
                        inf_engy.append(float(inf_temp.replace('+','E')))
                    elif('-' in inf_temp):
                        inf_engy.append(float(inf_temp.replace('-','E-')))

            col_excit.append(col_excit_row)
            #infinite energy

    ######################################################################
    #
    # end excitation and a values from the adf04, start looking for recomb
    #
    ######################################################################
        elif(tmp[0] == 'R' or tmp[0] == 'r'):
            tmp_inds = np.where(tmp!='')[0]
            recomb_transitions.append(np.array([int(tmp[tmp_inds[2]]),int(tmp[tmp_inds[1]])]))

            recomb_excit_row = np.zeros(len(adf04['input_file']['temp_grid']))
            for i in range(0,len(adf04['input_file']['temp_grid'])):
                if( '+' in tmp[tmp_inds[i+3]]):
                    recomb_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('+','E'))
                elif('-' in tmp[tmp_inds[i+3]]):
                    recomb_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('-','E-'))
            recomb_excit.append(recomb_excit_row)

    ######################################################################
    #
    # end recombination and values from the adf04, start looking for ionization
    #
    ######################################################################
        elif(tmp[0] == 'S' or tmp[0] =='s'):
            tmp_inds = np.where(tmp!='')[0]
            ion_transitions.append(np.array([int(tmp[tmp_inds[1]]),int(tmp[tmp_inds[2]])]))
            ion_excit_row = np.zeros(len(adf04['input_file']['temp_grid']))
            for i in range(0,len(adf04['input_file']['temp_grid'])):
                if( '+' in tmp[tmp_inds[i+3]]):
                    ion_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('+','E'))
                elif('-' in tmp[tmp_inds[i+3]]):
                    ion_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('-','E-'))
            ion_excit.append(ion_excit_row)
    adf04['rates']['excit'] = {}
    adf04['rates']['excit']['col_transitions'] = np.asarray(col_transitions)
    adf04['rates']['excit']['col_excit'] = np.asarray(col_excit)
    adf04['rates']['a_val'] = np.asarray(a_val)
    adf04['rates']['inf_engy'] = np.asarray(inf_engy)
    adf04['rates']['recomb'] = {}
    adf04['rates']['recomb']['recomb_transitions'] = np.asarray(recomb_transitions)
    adf04['rates']['recomb']['recomb_excit'] = np.asarray(recomb_excit)
    adf04['rates']['ioniz'] = {}
    adf04['rates']['ioniz']['ion_transitions'] = np.asarray(ion_transitions)
    adf04['rates']['ioniz']['ion_excit'] = np.asarray(ion_excit)

    adf04['input_file']['rates'] = {}
    adf04['input_file']['atomic'] = {}
    adf04['input_file']['rates'] = copy.deepcopy(adf04['rates'])
    adf04['input_file']['atomic'] = copy.deepcopy(adf04['atomic'])
    return adf04
