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
import sys
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
    f = open(file_name)
    ############################################################################
    #
    # reading the first line
    #
    ############################################################################
    try:
        first_line = f.readline()
        adf04['element'] = re.split('[+]',first_line)[0]
        first_linep = np.asarray(re.split('[ ]',re.split('[+]',first_line)[1]))
        first_linep_inds = np.where(first_linep !='')[0]
        adf04['charge_state'] = int(first_linep[first_linep_inds[0]])
        adf04['iz0'] = int(first_linep[first_linep_inds[1]])
        adf04['iz1'] = int(first_linep[first_linep_inds[2]])
        adf04['ion_pot'] = np.zeros(len(first_linep_inds)-3)
        adf04['ion_term'] = np.chararray(len(first_linep_inds)-3,itemsize=3)

        for i in range(0,len(adf04['ion_pot'])):
            adf04['ion_pot'][i] =float(
                re.split('[(]',first_linep[first_linep_inds[i+3]])[0])
            adf04['ion_term'][i] = re.split(
                '[)]',re.split('[(]',first_linep[first_linep_inds[i+3]])[1])[0]
    except:
        print("**Failed first line of the adf04**\n\n"+first_line+"\n\nfile is formatted incorrectly failed line.\nIt should be formatted as below (Neutral Tungsten example). \nIf more than one ionization metastable,\nadd more ionization potentials and terms after the first\nW+ 0        74         1       63427.7(6D3)")
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
                w.append(float(re.split('[)]',re.split('[(]',tmp[tmp_inds[config_stop+1]])[0])[0]))
                energy.append(float(tmp[tmp_inds[config_stop+2]]))
                zz = []
                zz1 = []
                for i in range(0,len(adf04['ion_term'])):
                    if(len(tmp_inds)>config_stop+3+i):
                        zz1.append(int(re.split('[}]',tmp[tmp_inds[config_stop+3+i]])[0][1]))
                        zz.append(float(re.split('[}]',tmp[tmp_inds[config_stop+3+i]])[1]))
                    else:
                        zz.append(-1.)
                        zz1.append(0)
                zpla.append(np.asarray(zz))
                zpla1.append(np.asarray(zz1))
        adf04['config'] = np.asarray(config)
        adf04['L'] = np.asarray(L)
        adf04['S'] = np.asarray(S)
        adf04['w'] = np.asarray(w)
        adf04['energy'] = np.asarray(energy)
        adf04['zpla'] = np.asarray(zpla)
        adf04['zpla1'] = np.asarray(zpla1)    
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
            elif('-' in temp_grid_tmp[i]):
                temp_tmp = temp_line_arr[temp_line_inds[i]].replace('-','E-')
            adf04_temp_grid.append(float(temp_tmp))
        adf04['temp_grid'] = np.asarray(adf04_temp_grid)
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
            col_excit_row = np.zeros(len(adf04['temp_grid']))
            for i in range(0,len(adf04['temp_grid'])):
                if( '+' in tmp[tmp_inds[i+3]]):
                    col_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('+','E'))
                elif('-' in tmp[tmp_inds[i+3]]):
                    col_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('-','E-'))
            col_excit.append(col_excit_row)
            #infinite energy
            if(len(tmp_inds)-3>=len(adf04['temp_grid'])):
                if( '+' in tmp[tmp_inds[-1]]):
                    inf_engy.append(float(tmp[tmp_inds[-1]].replace('+','E')))
                elif('-' in tmp[tmp_inds[-1]]):
                    inf_engy.append(float(tmp[tmp_inds[-1]].replace('-','E-')))
    ######################################################################
    #
    # end excitation and a values from the adf04, start looking for recomb
    #
    ######################################################################
        elif(tmp[0] == 'R' or tmp[0] == 'r'):
            tmp_inds = np.where(tmp!='')[0]
            recomb_transitions.append(np.array([int(tmp[tmp_inds[2]]),int(tmp[tmp_inds[1]])]))

            recomb_excit_row = np.zeros(len(adf04['temp_grid']))
            for i in range(0,len(adf04['temp_grid'])):
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
            ion_excit_row = np.zeros(len(adf04['temp_grid']))
            for i in range(0,len(adf04['temp_grid'])):
                if( '+' in tmp[tmp_inds[i+3]]):
                    ion_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('+','E'))
                elif('-' in tmp[tmp_inds[i+3]]):
                    ion_excit_row[i] = float(tmp[tmp_inds[i+3]].replace('-','E-'))
            ion_excit.append(ion_excit_row)

    adf04['col_transitions'] = np.asarray(col_transitions)
    adf04['col_excit'] = np.asarray(col_excit)
    adf04['a_val'] = np.asarray(a_val)
    adf04['inf_engy'] = np.asarray(inf_engy)
    adf04['recomb_transitions'] = np.asarray(recomb_transitions)
    adf04['recomb_excit'] = np.asarray(recomb_excit)
    adf04['ion_transitions'] = np.asarray(ion_transitions)
    adf04['ion_excit'] = np.asarray(ion_excit)

    return adf04
