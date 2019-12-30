import numpy as np
import re


fil = '../cbnm93#h_h0ls.dat'
f = open(fil)
adf17 = {}
adf17['input_file'] = {}

first_ln = f.readline()
first_ln_inds = [x.start() for x in re.finditer('=',first_ln)]

if(len(first_ln_inds) < 5):
    print('The first line of the file is in the wrong format')

##############################################################################################################
adf17['input_file']['sequence']    = re.search('\w+',first_ln[first_ln_inds[0]+1:]).group()
adf17['input_file']['nuc_chg']    = int(re.search('\w+',first_ln[first_ln_inds[1]+1:]).group())
adf17['input_file']['num_parent'] = int(re.search('\w+',first_ln[first_ln_inds[2]+1:]).group())
adf17['input_file']['num_dens']   = int(re.search('\w+',first_ln[first_ln_inds[3]+1:]).group())
adf17['input_file']['num_temp']   = int(re.search('\w+',first_ln[first_ln_inds[4]+1:]).group())


f.readline() #there is a blank line
second_ln = f.readline()

second_arr = re.split(' = ',second_ln)
adf17['input_file'][second_arr[0].replace(' ' ,'')] = ''

prev_dic_key = second_arr[0].replace(' ' ,'')
for i in range(1,len(second_arr)):
    tmp = re.split('  ', second_arr[i])
    adf17['input_file'][prev_dic_key] = int(tmp[0])
    if( i < len(second_arr) - 1):
        adf17['input_file'][tmp[1].replace(' ' ,'')] = ''        
        prev_dic_key = tmp[1].replace(' ' ,'')

#reading the densities
f.readline() #there is a blank line
dens_line = f.readline()
dens_arr = np.asarray([y.replace('D','e') for y in
            [x.strip(' ') for x in re.split('  ',re.split('=',dens_line)[1])[1:]]],dtype='float64')

while(len(dens_arr) < adf17['input_file']['num_dens']):
    dens_line = f.readline()
    dens_arr = np.append(dens_arr,np.asarray([y.replace('D','e') for y in
            [x.strip(' ') for x in re.split('  ',re.split('=',dens_line)[1])[1:]]],dtype='float64'))
adf17['input_file']['dens_grid'] = dens_arr
#reading the temperatures
temp_line = f.readline()
temp_arr = np.asarray([y.replace('D','e') for y in
            [x.strip(' ') for x in re.split('\d  ',re.split('=',temp_line)[1])]],dtype='float64')
while(len(temp_arr) < adf17['input_file']['num_temp']):
    temp_line = f.readline()
    temp_arr = np.append(temp_arr,np.asarray([y.replace('D','e') for y in
            [x.strip(' ') for x in re.split('\d  ',re.split('=',temp_line)[1])]],dtype='float64'))
adf17['input_file']['temp_grid'] = temp_arr

f.readline() #black line

parent_info_line = f.readline()
adf17['input_file']['num_parent']      = int(re.split('=',re.split('   ',parent_info_line)[0])[1])
adf17['input_file']['parent_term']     = re.split('=',re.split('   ',parent_info_line)[1])[1]
adf17['input_file']['num_spin_parent'] =  float(re.split('=',re.split('   ',parent_info_line)[2].replace('\n',''))[1])

f.readline()# ----------------line

spin_sys_line = f.readline()
spin_sys = float(re.split('\.', re.split('=',spin_sys_line.replace(' ', ''))[1])[0])
n_shells = int(re.split('=',spin_sys_line.replace(' ', ''))[2])

f.readline() # -------------------line




adf17['input_file']['spnsys_'+str(spin_sys)] = {}
adf17['input_file']['spnsys_'+str(spin_sys)]['excit'] = np.zeros((n_shells,n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['ion'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['rhs'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))

adf17['input_file']['spnsys_'+str(spin_sys)]['dir'] = {}
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['excit'] = np.zeros((n_shells,n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['ion'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_three_body'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_diele'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_rad'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))
adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['cxe'] = np.zeros((n_shells,len(adf17['input_file']['temp_grid']),len(adf17['input_file']['dens_grid'])))





for d in range(0,len(adf17['input_file']['dens_grid'])):
    for t in range(0,len(adf17['input_file']['temp_grid'])):

        ttmp = f.readline()
        s = re.split('  ',ttmp)[1].replace('D','e')
        adf17['input_file']['spnsys_'+str(spin_sys)]['excit'][0,:,t,d] = np.asarray([s[i:i+17] for i in range(0, len(s), 17)],dtype='float64')
        adf17['input_file']['spnsys_'+str(spin_sys)]['ion'][0,t,d] = float(re.split('  ',ttmp)[2].replace('D','e'))
        adf17['input_file']['spnsys_'+str(spin_sys)]['rhs'][0,t,d] = float(re.split('  ',ttmp)[-1].replace('\n','').replace('D','e'))


        for i in range(1,n_shells):
            ttmp = f.readline()
            s = re.split('  ',ttmp)[7].replace('D','e')
            adf17['input_file']['spnsys_'+str(spin_sys)]['excit'][i,:,t,d] = np.asarray([s[i:i+17] for i in range(0, len(s), 17)],dtype='float64')
            adf17['input_file']['spnsys_'+str(spin_sys)]['ion'][i,t,d] = float(re.split('  ',ttmp)[8].replace('D','e'))
            adf17['input_file']['spnsys_'+str(spin_sys)]['rhs'][i,t,d] = float(re.split('  ',ttmp)[-1].replace('\n','').replace('D','e'))

        ttmp = f.readline()
        s = re.split('  ',ttmp)[4].replace('D','e')
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['excit'][0,:,t,d] = np.asarray([s[i:i+17] for i in range(0, len(s), 17)],dtype='float64')
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['ion'][0,t,d] = float(re.split('  ',ttmp)[5].replace('D','e'))
        ss = re.split('  ',ttmp)[-1].replace('D','e').replace('\n','')
        ss_arr = np.asarray([ss[i:i+17] for i in range(0, len(ss), 17)],dtype='float64')
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_three_body'][0,t,d] =ss_arr[0]
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_diele'][0,t,d] = ss_arr[1]
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_rad'][0,t,d] = ss_arr[2]
        adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['cxe'][0,t,d] = ss_arr[3]


        for i in range(1,n_shells):
            ttmp = f.readline()
            s = re.split('  ',ttmp)[7].replace('D','e')
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['excit'][i,:,t,d] = np.asarray([s[i:i+17] for i in range(0, len(s), 17)],dtype='float64')
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['ion'][i,t,d] = float(re.split('  ',ttmp)[8].replace('D','e'))
            ss = re.split('  ',ttmp)[-1].replace('D','e').replace('\n','')
            ss_arr = np.asarray([ss[i:i+17] for i in range(0, len(ss), 17)],dtype='float64')
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_three_body'][i,t,d] =ss_arr[0]
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_diele'][i,t,d] = ss_arr[1]
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['recomb_rad'][i,t,d] = ss_arr[2]
            adf17['input_file']['spnsys_'+str(spin_sys)]['dir']['cxe'][i,t,d] = ss_arr[3]

        f.readline()#pbr
        f.readline()#blank line
