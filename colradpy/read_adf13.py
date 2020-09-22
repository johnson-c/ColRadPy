import numpy as np
import re



#fil = '/home/curtis/adas/adas/adas/adf13/sxb96#c/sxb96#c_vsu#c1.dat'




def read_adf13(fil):
    f  = open(fil)

    adf13 = {}
    adf13['input_file'] = {}



    #reading the first line
    q = f.readline()
    adf13['input_file']['num_trans'] = int(re.split('/',q)[0])
    adf13['input_file']['element'] = re.split( '\+',re.split('IONI',re.split('/',q)[1])[0])[0]
    adf13['input_file']['charge_state'] = int(re.split( '\+',re.split('IONI',re.split('/',q)[1])[0])[1])



    #reading the transition line
    grid_line = f.readline()



    tmp = 'A' in re.split('\s+',grid_line)[1]

    if(tmp):
        adf13['input_file']['num_temp'] = int(re.split('\s+',grid_line)[3])
        adf13['input_file']['num_dens'] = int(re.split('\s+',grid_line)[2])
        adf13['input_file']['orig_file'] = re.split('\s+',grid_line)[6]
        adf13['input_file']['meta_indx'] = re.split('\s+',grid_line)[12]

    else:
        adf13['input_file']['num_temp'] = int(re.split('\s+',grid_line)[4])
        adf13['input_file']['num_dens'] = int(re.split('\s+',grid_line)[3])
        adf13['input_file']['orig_file'] = re.split('\s+',grid_line)[7]
        adf13['input_file']['meta_indx'] = re.split('\s+',grid_line)[13]






    adf13['input_file']['wave_air'] = []
    adf13['input_file']['sxb'] = []

    while 'C-' not in grid_line:
        if(tmp):

            adf13['input_file']['wave_air'].append(float(re.split('A',re.split('\s+',grid_line)[1])[0])) 
        else:
            adf13['input_file']['wave_air'].append(float(re.split('\s+',grid_line)[1]))
        #read in the density grid
        adf13['input_file']['dens_grid'] = np.array([])
        while adf13['input_file']['dens_grid'].size <adf13['input_file']['num_dens']:
            den_line = f.readline()
            den_line = den_line.rstrip()
            adf13['input_file']['dens_grid'] = np.append(adf13['input_file']['dens_grid'], np.array(list(map(float,re.split('\s+',den_line)[1:])))  )


        #read in the temperature grid
        adf13['input_file']['temp_grid'] = np.array([])
        while adf13['input_file']['temp_grid'].size <adf13['input_file']['num_temp']:

            temp_line = f.readline()
            temp_line = temp_line.rstrip()
            adf13['input_file']['temp_grid'] = np.append(adf13['input_file']['temp_grid'], np.array(list(map(float,re.split('\s+',temp_line)[1:])))  )

        grid_line = f.readline()
        t = []
        while 'A' not in grid_line and 'C-' not in grid_line:

            t.append(np.array(list(map(float,re.split('\s+',grid_line.rstrip())[1:]))))
            grid_line = f.readline()
            
        adf13['input_file']['sxb'].append(np.reshape(np.concatenate(t),(adf13['input_file']['num_dens'],adf13['input_file']['num_temp'])))

    adf13['input_file']['wave_air'] = np.asarray(adf13['input_file']['wave_air'])
    adf13['input_file']['sxb'] = np.asarray(adf13['input_file']['sxb'])

    adf13['input_file']['sxb'] = adf13['input_file']['sxb'].transpose(0,2,1)
                                          
    return adf13
