import numpy as np
import os
from nist_read_txt import *
import requests
import re



auto_get_fe = True
auto_get_kr = False
auto_get_mo = False


def get_from_open_adas(t, path='../atomic/adf04_downloaded_from_open_adas/',add_ionization_blind = False,element = '', charge_state = ''):
    
    a = requests.get(t)
    file_name = re.split('\.',re.split('=',a.headers['Content-Disposition'])[1])[0] + '.mod'
    ttmp = re.sub('\s+\n', '\n',a.text)#remove trailing white space if its in the file
    if(add_ionization_blind):
        ttmp = re.sub('{X}', '{1}1.00',ttmp )


    if( '(' not in re.split('\n',ttmp[0:200])[0]):#check if the people even added the term
        ttmp = re.sub(r'\n',r'(   )\n',ttmp,count=1)        
    if(charge_state != ''):
        ter = get_nist_txt('/home/curtis/git/high-z_line_id/atomic/nist_energies/',element,charge_state)[0]['term']
        ttmp = re.sub(r'\(\s+\)','(' + ter + ')',ttmp)
    else:
        ttmp = re.sub(r'\(\s+\)','(1S)',ttmp)        
        
    open(path+file_name,'w').write(ttmp)



url = 'http://open.adas.ac.uk/download/adf04/'

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



#################################
#
# Iron data set start
#
#################################

fe_names_adas  = np.array([
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


fe_names_local = np.array(['kr/ssh41_cs_ic#fe0.mod',       'kr/ssh41_cs_ic#fe1.mod',
                           'kr/crlike_nrb13#fe2.mod',      'kr/ssh41_cs_ic#fe3.mod',
                           'kr/ssh41_cs_ic#fe4.mod',       'kr/ssh41_cs_ic#fe5.mod',
                           'kr/calike_mcw07#fe6.mod',      'klike_gdz14#fe7.mod',
                           'kr/arlike_dgz14#fe8.mod',      'kr/ssh41_cs_ic#fe9.mod',
                           'kr/ssh41_cs_ic#fe10.mod',      'kr/ssh41_cs_ic#fe11.mod',
                           'kr/ssh41_cs_ic#fe12.mod',      'kr/ssh41_cs_ic#fe13.mod',
                           'kr/ssh41_cs_ic#fe14.mod',      'kr/ssh41_cs_ic#fe15.mod',
                           'kr/ssh41_cs_ic#fe16.mod',      'kr/ssh41_cs_ic#fe17.mod',
                           'kr/ssh41_cs_ic#fe18.mod',      'kr/ssh41_cs_ic#fe19.mod',
                           'kr/ssh41_cs_ic#fe20.mod',      'kr/ssh41_cs_ic#fe21.mod',
                           'kr/ssh41_cs_ic#fe22.mod',      'kr/ssh41_cs_ic#fe23.mod',
                           'kr/ssh41_cs_ic#fe24.mod',      'kr/ssh41_cs_ic#fe25.mod'])



mo_names_local = np.array(['mo/ic#mo0.mod',            'mo/ic#mo1.mod',
                           'mo/ic#mo2.mod',            'mo/ic#mo3.mod',
                           'mo/ic#mo4.mod',            'mo/ic#mo5.mod',
                           'mo/ic#mo6.mod',            'mo/ic#mo7.mod',
                           'mo/ic#mo8.mod',            'mo/ic#mo9.mod',
                           'mo/ic#mo10.mod',           'mo/ic#mo11.mod',
                           'mo/ic#mo12.mod',           'mo/ic#mo13.mod',
                           'mo/ic#mo14.mod',           'mo/ic#mo15.mod',
                           'mo/ic#mo16.mod',           'mo/ic#mo17.mod',
                           'mo/ic#mo18.mod',           'mo/ic#mo19.mod',
                           'mo/ic#mo20.mod',           'mo/ic#mo21.mod',
                           'mo/ic#mo22.mod',           'mo/ic#mo23.mod',
                           'mo/ic#mo24.mod',           'mo/ic#mo25.mod',
                           'mo/ic#mo26.mod',           'mo/ic#mo27.mod',
                           'mo/ic#mo28.mod',           'mo/ic#mo29.mod',
                           'mo/ic#mo30.mod',           'mo/ic#mo31.mod',
                           'mo/ic#mo32.mod',           'mo/copsm#f_sm#mo33.mod',
                           'mo/ic#mo34.mod',           'mo/ic#mo35.mod',
                           'mo/ic#mo36.mod',           'mo/ic#mo37.mod',
                           'mo/copsm#be_sm#mo38.mod',  'mo/ic#mo39.mod',
                           'mo/ic#mo40.mod',           'mo/ic#mo41.mod'])

kr_names_local = np.array(['kr/ssh41_cs_ic#kr0.mod',       'kr/ssh41_cs_ic#kr1.mod',
                           'kr/ssh41_cs_ic#kr2.mod',       'kr/ssh41_cs_ic#kr3.mod',
                           'kr/ssh41_cs_ic#kr4.mod',       'kr/ssh41_cs_ic#kr5.mod',
                           'kr/ssh41_cs_ic#kr6.mod',       'kr/ssh41_cs_ic#kr7.mod',
                           'kr/ssh41_cs_ic#kr8.mod',       'kr/ssh41_cs_ic#kr9.mod',
                           'kr/ssh41_cs_ic#kr10.mod',      'kr/ssh41_cs_ic#kr11.mod',
                           'kr/ssh41_cs_ic#kr12.mod',      'kr/ssh41_cs_ic#kr13.mod',
                           'kr/ssh41_cs_ic#kr14.mod',      'kr/ssh41_cs_ic#kr15.mod',
                           'kr/ssh41_cs_ic#kr16.mod',      'kr/ssh41_cs_ic#kr17.mod',
                           'kr/ssh41_cs_ic#kr18.mod',      'kr/ssh41_cs_ic#kr19.mod',
                           'kr/ssh41_cs_ic#kr20.mod',      'kr/ssh41_cs_ic#kr21.mod',
                           'kr/ssh41_cs_ic#kr22.mod',      'kr/ssh41_cs_ic#kr23.mod',
                           'kr/mglike_lfm14#kr24.mod',     'kr/nalike_lgy09#kr25.mod',
                           'kr/nelike_dcg08#kr26.mod',     'kr/flike_mcw06#kr27.mod',
                           'kr/ssh41_cs_ic#kr28.mod',      'kr/ssh41_cs_ic#kr29.mod',
                           'kr/clike_jm19#kr30.mod',       'kr/blike_lgy12#kr31.mod',
                           'kr/belike_lfm14#kr32.mod',     'kr/lilike_lgy10#kr33.mod',
                           'kr/helike_adw05#kr34.mod',     'kr/ssh41_cs_ic#kr35.mod'])











if(auto_get_fe):
    fe_path = '../atomic/adf04_downloaded_from_open_adas/fe/'
    if not os.path.exists(fe_path):
        os.makedirs(fe_path)
    for i in range(0,len(fe_names_adas)):
        if( i < len(fe_names_adas) -1):
            get_from_open_adas(t = url+fe_names_adas[i],path = fe_path, add_ionization_blind=True,element='Fe',charge_state=i+2)
        else:
            get_from_open_adas(t = url+fe_names_adas[i],path = fe_path, add_ionization_blind=True)
    

#################################
#
# Iron data set end
#
#################################









    
    
kr_names_adas = np.array(['copssh][36/ssh41_cs_ic][kr0.dat',          'copssh][36/ssh41_cs_ic][kr1.dat',
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






if(auto_get_kr):
    kr_path = '../atomic/adf04_downloaded_from_open_adas/kr/'
    if not os.path.exists(kr_path):
        os.makedirs(kr_path)
    for i in range(0,len(kr_names_adas)):
        if( i < len(kr_names_adas) -1):
            get_from_open_adas(t = url+kr_names_adas[i],path = kr_path, add_ionization_blind=True,element='Kr',charge_state=i+2)
        else:
            get_from_open_adas(t = url+kr_names_adas[i],path = kr_path, add_ionization_blind=True)
    











mo_names_adas  = np.array([
          'copmm][42/ic][mo0.dat',               'copmm][42/ic][mo1.dat',
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




if(auto_get_mo):
    mo_path = '../atomic/adf04_downloaded_from_open_adas/mo/'
    if not os.path.exists(mo_path):
        os.makedirs(mo_path)
    for i in range(0,3):#2,len(mo_names_adas)):
        if( i < len(mo_names_adas) -1):
            get_from_open_adas(t = url+mo_names_adas[i],path = mo_path, add_ionization_blind=True,element='Mo',charge_state=i+2)
        else:
            get_from_open_adas(t = url+mo_names_adas[i],path = mo_path, add_ionization_blind=True)
                                                   

                                                   



                                                                         
