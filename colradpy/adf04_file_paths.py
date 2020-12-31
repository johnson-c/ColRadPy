import numpy as np
import os
from colradpy.nist_read_txt import *
import requests
import re


def get_from_open_adas(tt,charge_state, path='../atomic/adf04_downloaded_from_open_adas/',add_ionization_blind = True):

    url = url = 'http://open.adas.ac.uk/download/adf04/'

    add_ionization_blind = True
    adas_paths = {}
    adas_paths['Fe']  = np.array([
        'copssh][26/ssh41_cs_ic][fe0.dat',
        'copssh][26/ssh41_cs_ic][fe1.dat',
        'crlike/crlike_nrb13][fe2.dat',
        'copssh][26/ssh41_cs_ic][fe3.dat',
        'copssh][26/ssh41_cs_ic][fe4.dat',
        'copssh][26/ssh41_cs_ic][fe5.dat',
        'calike/calike_mcw07][fe6.dat',
        'klike/klike_gdz14][fe7.dat',
        'arlike/arlike_dgz14][fe8.dat',
        'copssh][26/ssh41_cs_ic][fe9.dat',
        'copssh][26/ssh41_cs_ic][fe10.dat',
        'copssh][26/ssh41_cs_ic][fe11.dat',
        'copssh][26/ssh41_cs_ic][fe12.dat',
        'copssh][26/ssh41_cs_ic][fe13.dat',
        'copssh][26/ssh41_cs_ic][fe14.dat',
        'copssh][26/ssh41_cs_ic][fe15.dat',
        'copssh][26/ssh41_cs_ic][fe16.dat',
        'copssh][26/ssh41_cs_ic][fe17.dat',
        'copssh][26/ssh41_cs_ic][fe18.dat',
        'copssh][26/ssh41_cs_ic][fe19.dat',
        'copssh][26/ssh41_cs_ic][fe20.dat',
        'copssh][26/ssh41_cs_ic][fe21.dat',
        'copssh][26/ssh41_cs_ic][fe22.dat',
        'copssh][26/ssh41_cs_ic][fe23.dat',
        'copssh][26/ssh41_cs_ic][fe24.dat',    
        'copssh][26/ssh41_cs_ic][fe25.dat'])


    adas_paths['Kr'] = np.array(['copssh][36/ssh41_cs_ic][kr0.dat',          'copssh][36/ssh41_cs_ic][kr1.dat',
                          'copssh][36/ssh41_cs_ic][kr2.dat',          'copssh][36/ssh41_cs_ic][kr3.dat',
                          'copssh][36/ssh41_cs_ic][kr4.dat',          'copssh][36/ssh41_cs_ic][kr5.dat',
                          'copssh][36/ssh41_cs_ic][kr6.dat',          'copssh][36/ssh41_cs_ic][kr7.dat',
                          'copssh][36/ssh41_cs_ic][kr8.dat',          'copssh][36/ssh41_cs_ic][kr9.dat',
                          'copssh][36/ssh41_cs_ic][kr10.dat',         'copssh][36/ssh41_cs_ic][kr11.dat',
                          'copssh][36/ssh41_cs_ic][kr12.dat',         'copssh][36/ssh41_cs_ic][kr13.dat',
                          'copssh][36/ssh41_cs_ic][kr14.dat',         'copssh][36/ssh41_cs_ic][kr15.dat',
                          'copssh][36/ssh41_cs_ic][kr16.dat',         'copssh][36/ssh41_cs_ic][kr17.dat',
                          'copssh][36/ssh41_cs_ic][kr18.dat',         'copssh][36/ssh41_cs_ic][kr19.dat',
                          'copssh][36/ssh41_cs_ic][kr20.dat',         'copssh][36/ssh41_cs_ic][kr21.dat',
                          'copssh][36/ssh41_cs_ic][kr22.dat',         'copssh][36/ssh41_cs_ic][kr23.dat',
                          'copaw][mg/mglike_lfm14][kr24.dat',         'copaw][na/nalike_lgy09][kr25.dat',
                          'nelike/nelike_dcg08][kr26.dat',            'copaw][f/flike_mcw06][kr27.dat',
                          'copssh][36/ssh41_cs_ic][kr28.dat',         'copssh][36/ssh41_cs_ic][kr29.dat',
                          'copaw][c/clike_jm19][kr30.dat',            'copaw][b/blike_lgy12][kr31.dat',
                          'copaw][be/belike_lfm14][kr32.dat',         'copaw][li/lilike_lgy10][kr33.dat',
                          'copaw][he/helike_adw05][kr34.dat',         'copssh][36/ssh41_cs_ic][kr35.dat'])



    adas_paths['Mo']= np.array(['copmm][42/ic][mo0.dat',               'copmm][42/ic][mo1.dat',
                                'copmm][42/ic][mo2.dat',               'copmm][42/ic][mo3.dat',
                                'copmm][42/ic][mo4.dat',               'copmm][42/ic][mo5.dat',
                                'copmm][42/ic][mo6.dat',               'copmm][42/ic][mo7.dat',
                                'copmm][42/ic][mo8.dat',               'copmm][42/ic][mo9.dat',
                                'copmm][42/ic][mo10.dat',              'copmm][42/ic][mo11.dat',
                                'copmm][42/ic][mo12.dat',              'copmm][42/ic][mo13.dat',
                                'copmm][42/ic][mo14.dat',              'copmm][42/ic][mo15.dat',
                                'copmm][42/ic][mo16.dat',              'copmm][42/ic][mo17.dat',
                                'copmm][42/ic][mo18.dat',              'copmm][42/ic][mo19.dat',
                                'copmm][42/ic][mo20.dat',              'copmm][42/ic][mo21.dat',
                                'copmm][42/ic][mo22.dat',              'copmm][42/ic][mo23.dat',
                                'copmm][42/ic][mo24.dat',              'copmm][42/ic][mo25.dat',
                                'copmm][42/ic][mo26.dat',              'copmm][42/ic][mo27.dat',
                                'copmm][42/ic][mo28.dat',              'copmm][42/ic][mo29.dat',
                                'copmm][42/ic][mo30.dat',              'copmm][42/ic][mo31.dat',
                                'copmm][42/ic][mo32.dat',              'copsm][f/copsm][f_sm][mo33.dat',
                                'copmm][42/ic][mo34.dat',              'copmm][42/ic][mo35.dat',
                                'copmm][42/ic][mo36.dat',              'copmm][42/ic][mo37.dat',
                                'copsm][be/copsm][be_sm][mo38.dat',    'copmm][42/ic][mo39.dat',
                                'copmm][42/ic][mo40.dat',              'copmm][42/ic][mo41.dat'])





    t = url + adas_paths[tt][charge_state]
    
    
    a = requests.get(t)
    file_name = re.split('\.',re.split('=',a.headers['Content-Disposition'])[1])[0] + '.mod'
    ttmp = re.sub('\s+\n', '\n',a.text)#remove trailing white space if its in the file
    if(add_ionization_blind):
        ttmp = re.sub('{X}', '{1}1.00',ttmp )


    if( '(' not in re.split('\n',ttmp[0:200])[0]):#check if the people even added the term
        ttmp = re.sub(r'\n',r'(   )\n',ttmp,count=1)
        
    if(charge_state < len(adas_paths[tt]) -1):

        base_path = os.path.dirname(os.path.realpath(__file__))
        base_path = base_path[0:len(base_path) - 8] + 'atomic/'

        ter = get_nist_txt(base_path+'nist_energies/',tt.lower(),charge_state+2)[0]['term']
        ttmp = re.sub(r'\(\s+\)','(' + ter + ')',ttmp)
    else:
        ttmp = re.sub(r'\(\s+\)','(1S)',ttmp)


    local_dict = return_local_adf04_dict()
    

    if not os.path.exists(base_path+'adf04_downloaded_from_open_adas/'+tt):
        os.makedirs(base_path+'adf04_downloaded_from_open_adas/'+tt)

        
    open(local_dict[tt][charge_state],'w').write(ttmp)



def return_local_adf04_dict():
    
    paths = {}

    base_path = os.path.dirname(os.path.realpath(__file__))
    base_path = base_path[0:len(base_path) - 8] + 'atomic/adf04_downloaded_from_open_adas/'
        
    paths['Fe'] =  np.array([base_path + 'Fe/ssh41_cs_ic#fe0.mod',       base_path + 'Fe/ssh41_cs_ic#fe1.mod',
                             base_path + 'Fe/crlike_nrb13#fe2.mod',      base_path + 'Fe/ssh41_cs_ic#fe3.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe4.mod',       base_path + 'Fe/ssh41_cs_ic#fe5.mod',
                             base_path + 'Fe/calike_mcw07#fe6.mod',      'klike_gdz14#fe7.mod',
                             base_path + 'Fe/arlike_dgz14#fe8.mod',      base_path + 'Fe/ssh41_cs_ic#fe9.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe10.mod',      base_path + 'Fe/ssh41_cs_ic#fe11.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe12.mod',      base_path + 'Fe/ssh41_cs_ic#fe13.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe14.mod',      base_path + 'Fe/ssh41_cs_ic#fe15.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe16.mod',      base_path + 'Fe/ssh41_cs_ic#fe17.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe18.mod',      base_path + 'Fe/ssh41_cs_ic#fe19.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe20.mod',      base_path + 'Fe/ssh41_cs_ic#fe21.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe22.mod',      base_path + 'Fe/ssh41_cs_ic#fe23.mod',
                             base_path + 'Fe/ssh41_cs_ic#fe24.mod',      base_path + 'Fe/ssh41_cs_ic#fe25.mod'])


    
    paths['Mo'] = np.array(['Mo/ic#mo0.mod',            base_path + 'Mo/ic#mo1.mod',
                           base_path + 'Mo/ic#mo2.mod',            base_path + 'Mo/ic#mo3.mod',
                           base_path + 'Mo/ic#mo4.mod',            base_path + 'Mo/ic#mo5.mod',
                           base_path + 'Mo/ic#mo6.mod',            base_path + 'Mo/ic#mo7.mod',
                           base_path + 'Mo/ic#mo8.mod',            base_path + 'Mo/ic#mo9.mod',
                           base_path + 'Mo/ic#mo10.mod',           base_path + 'Mo/ic#mo11.mod',
                           base_path + 'Mo/ic#mo12.mod',           base_path + 'Mo/ic#mo13.mod',
                           base_path + 'Mo/ic#mo14.mod',           base_path + 'Mo/ic#mo15.mod',
                           base_path + 'Mo/ic#mo16.mod',           base_path + 'Mo/ic#mo17.mod',
                           base_path + 'Mo/ic#mo18.mod',           base_path + 'Mo/ic#mo19.mod',
                           base_path + 'Mo/ic#mo20.mod',           base_path + 'Mo/ic#mo21.mod',
                           base_path + 'Mo/ic#mo22.mod',           base_path + 'Mo/ic#mo23.mod',
                           base_path + 'Mo/ic#mo24.mod',           base_path + 'Mo/ic#mo25.mod',
                           base_path + 'Mo/ic#mo26.mod',           base_path + 'Mo/ic#mo27.mod',
                           base_path + 'Mo/ic#mo28.mod',           base_path + 'Mo/ic#mo29.mod',
                           base_path + 'Mo/ic#mo30.mod',           base_path + 'Mo/ic#mo31.mod',
                           base_path + 'Mo/ic#mo32.mod',           base_path + 'Mo/copsm#f_sm#mo33.mod',
                           base_path + 'Mo/ic#mo34.mod',           base_path + 'Mo/ic#mo35.mod',
                           base_path + 'Mo/ic#mo36.mod',           base_path + 'Mo/ic#mo37.mod',
                           base_path + 'Mo/copsm#be_sm#mo38.mod',  base_path + 'Mo/ic#mo39.mod',
                           base_path + 'Mo/ic#mo40.mod',           base_path + 'Mo/ic#mo41.mod'])

    
    paths['Kr'] =  np.array([base_path + 'Kr/ssh41_cs_ic#kr0.mod',       base_path + 'Kr/ssh41_cs_ic#kr1.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr2.mod',       base_path + 'Kr/ssh41_cs_ic#kr3.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr4.mod',       base_path + 'Kr/ssh41_cs_ic#kr5.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr6.mod',       base_path + 'Kr/ssh41_cs_ic#kr7.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr8.mod',       base_path + 'Kr/ssh41_cs_ic#kr9.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr10.mod',      base_path + 'Kr/ssh41_cs_ic#kr11.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr12.mod',      base_path + 'Kr/ssh41_cs_ic#kr13.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr14.mod',      base_path + 'Kr/ssh41_cs_ic#kr15.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr16.mod',      base_path + 'Kr/ssh41_cs_ic#kr17.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr18.mod',      base_path + 'Kr/ssh41_cs_ic#kr19.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr20.mod',      base_path + 'Kr/ssh41_cs_ic#kr21.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr22.mod',      base_path + 'Kr/ssh41_cs_ic#kr23.mod',
                           base_path + 'Kr/mglike_lfm14#kr24.mod',     base_path + 'Kr/nalike_lgy09#kr25.mod',
                           base_path + 'Kr/nelike_dcg08#kr26.mod',     base_path + 'Kr/flike_mcw06#kr27.mod',
                           base_path + 'Kr/ssh41_cs_ic#kr28.mod',      base_path + 'Kr/ssh41_cs_ic#kr29.mod',
                           base_path + 'Kr/clike_jm19#kr30.mod',       base_path + 'Kr/blike_lgy12#kr31.mod',
                           base_path + 'Kr/belike_lfm14#kr32.mod',     base_path + 'Kr/lilike_lgy10#kr33.mod',
                           base_path + 'Kr/helike_adw05#kr34.mod',     base_path + 'Kr/ssh41_cs_ic#kr35.mod'])


    return paths





'''
h_names_adas = np.array(['adas][1/ha00_ls][h0.dat'])
h_names_local = np.array(['ha00_ls#h0.mod'])


he_arr = np.array(['adas][2/mom97_ls][he0.dat','adas][2/mom97_ls][he1.dat'])
he_file = np.array(['he0_adf04','he1_adf04'])

li_arr = np.array(['adas][3/cpb02_ls][li0.dat','adas][3/cpb02_ls][li1.dat','adas][3/cpb02_ls][li2.dat'])
li_file = np.array(['li0_adf04','li1_adf04','li2_adf04'])

be_arr = np.array(['adas][4/cpb03_ls][be0.dat','adas][4/cpb03_ls][be1.dat','adas][4/cpb03_ls][be2.dat','adas][4/cpb03_ls][be3.dat'])
be_file = np.array(['be0_adf04','be1_adf04','be2_adf04','be3_adf04'])
                  

be_arr = np.array(['adas][4/cpb03_ls][be0.dat','adas][4/cpb03_ls][be1.dat','adas][4/cpb03_ls][be2.dat','adas][4/cpb03_ls][be3.dat'])
be_file = np.array(['be0_adf04','be1_adf04','be2_adf04','be3_adf04'])

b_arr = np.array(['cophps][b/dw/ls][b0.dat','copaw][be/belike_lfm14][b1.dat','copaw][li/lilike_lgy10][b2.dat','cophps][he/dw/ls][b3.dat','hlike/hlike_cpb02][b4.dat'])
b_file = np.array(['b0_adf04','b1_adf04','b2_adf04','b3_adf04'])


c_arr  = np.array(['adas][6/mom97_ls][c0.dat','adas][6/mom97_ls][c1.dat','adas][6/mom97_ls][c2.dat','adas][6/mom97_ls][c3.dat','adas][6/mom97_ls][c4.dat','adas][6/mom97_n][c5.dat'])
c_file = np.array(['c0_adf04','c1_adf04','c2_adf04','c3_adf04','c4_adf04','c5_adf04'])
'''
                                                                         
