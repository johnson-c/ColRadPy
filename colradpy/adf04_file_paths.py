import urllib.request
import numpy as np
import os

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

https://open.adas.ac.uk/detail/adf04/

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


fe_names_local = np.array(['ssh41_cs_ic#fe0.mod',
                    'ssh41_cs_ic#fe1.mod',
                    'crlike_nrb13#fe2.mod',
                    'ssh41_cs_ic#fe3.mod',
                    'ssh41_cs_ic#fe4.mod',
                    'ssh41_cs_ic#fe5.mod',
                    'calike_mcw07#fe6.mod',
                    'klike_gdz14#fe7.mod',
                    'arlike_dgz14#fe8.mod',
                    'ssh41_cs_ic#fe9.mod',
                    'ssh41_cs_ic#fe10.mod',
                    'ssh41_cs_ic#fe11.mod',
                    'ssh41_cs_ic#fe12.mod',
                    'ssh41_cs_ic#fe13.mod',
                    'ssh41_cs_ic#fe14.mod',
                    'ssh41_cs_ic#fe15.mod',
                    'ssh41_cs_ic#fe16.mod',
                    'ssh41_cs_ic#fe17.mod',
                    'ssh41_cs_ic#fe18.mod',
                    'ssh41_cs_ic#fe19.mod',
                    'ssh41_cs_ic#fe20.mod',
                    'ssh41_cs_ic#fe21.mod',
                    'ssh41_cs_ic#fe22.mod',
                    'ssh41_cs_ic#fe23.mod',
                    'ssh41_cs_ic#fe24.mod',
                    'ssh41_cs_ic#fe25.mod'])

'''

iron_path = '../atomic/adf04_downloaded_from_open_adas/fe'
if not os.path.exists(iron_path):
    os.makedirs(iron_path)

for i in range(0,len(fe_arr)):
    get_from_open_adas(t = url+fe_arr[i],path = iron_path add_ionization_blind=True)
'''


#################################
#
# Iron data set end
#
#################################




















def get_from_open_adas(t, path='../atomic/adf04_downloaded_from_open_adas/',add_ionization_blind = False,):
    
    a = requests.get(t)
    file_name = re.split('\.',re.split('=',a.headers['Content-Disposition'])[1])[0] + '.mod'
    ttmp = re.sub('\s+\n', '\n',a.text)#remove trailing white space if its in the file
    if(add_ionization_blind):
        ttmp = re.sub('{X}', '{1}1.00',ttmp )
        
    open(path+file_name,'w').write(ttmp)



    

