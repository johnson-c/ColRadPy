import numpy as np
import re

f = open('/home/curtis/Downloads/scd96r_he.dat')
adf11 = {}
adf11['input_file'] = {}



#reading the first line
tmp = re.findall('(\d+)',f.readline())
adf11['input_file']['nuc_charge'] = int(tmp[0])
adf11['input_file']['num_dens'] = int(tmp[1])
adf11['input_file']['num_temp'] = int(tmp[2])
adf11['input_file']['charge_min'] = int(tmp[3])
adf11['input_file']['charge_max'] = int(tmp[4])

f.readline() #reading '-------------'
adf11['input_file']['metas'] = np.array(list(   #metastables
               map(int,re.findall('(\d+)',f.readline()))))
f.readline() #reading '---------------'

#read in the density grid
adf11['input_file']['dens'] = np.array([])
while adf11['input_file']['dens'].size <adf11['input_file']['num_dens']:
    adf11['input_file']['dens'] = np.append(adf11['input_file']['dens'],np.array(list(map(float,re.findall('(.\d*\.\d+)',f.readline())))))
#read in the temperature grid
adf11['input_file']['temp'] = np.array([])
while adf11['input_file']['temp'].size <adf11['input_file']['num_temp']:
    adf11['input_file']['temp'] = np.append(adf11['input_file']['temp'],np.array(list(map(float,re.findall('(.\d*\.\d+)',f.readline())))))
    
