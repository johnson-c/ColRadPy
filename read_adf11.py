import numpy as np
import re


def read_adf11(fil):
    """ Reads data formatted in adf11 format. This format is generally used for GCR coefficients.
        The data is stored in arrays of log_10.
        Creates a dictionary to hold the data
        Args:
          :param fil: The file path to the input file
      :type fil: string
    """
    
    f = open(fil)
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

    #setting up the GCR stages
    if(adf11['input_file']['nuc_charge'] == len(adf11['input_file']['metas'])-1):
        num_stages = adf11['input_file']['nuc_charge']
    else:
        num_stages = len(adf11['input_file']['metas'])
    for i in range(0,num_stages):
        adf11['input_file'][str(i)] = {}
        if( 'scd' in fil or 'acd' in fil):
            adf11['input_file'][str(i)] = np.zeros((adf11['input_file']['metas'][i],
                                                     adf11['input_file']['metas'][i+1],
                                             len(adf11['input_file']['temp']),
                                             len(adf11['input_file']['dens'])))
        if( 'qcd' in fil ):
            adf11['input_file'][str(i)] = np.zeros((adf11['input_file']['metas'][i],
                                                     adf11['input_file']['metas'][i],
                                             len(adf11['input_file']['temp']),
                                             len(adf11['input_file']['dens'])))
        if('xcd' in fil):
            adf11['input_file'][str(i)] = np.zeros((adf11['input_file']['metas'][i+1],
                                                     adf11['input_file']['metas'][i+1],
                                             len(adf11['input_file']['temp']),
                                             len(adf11['input_file']['dens'])))
    #Reading the GCR value portion
    gcr_line = f.readline()
    ii = 0
    while 'C' not in gcr_line:
        #look for the stage identifying line
        if('---' in gcr_line):
            dens_count = 0
            temp_count = 0
            stage_id = np.array(list(map(int,re.findall('( \d+ )',gcr_line))))

            if('scd' in fil or 'acd' in fil):
                tmp = stage_id[0]
                stage_id[0] = stage_id[1]
                stage_id[1] = tmp

            ii = ii+1
        else:
            gcr_vals = np.array(list(map(float,re.findall('(-\d+.\d+)',gcr_line))))
            adf11['input_file'][str(stage_id[2]-1)][stage_id[0]-1,stage_id[1]-1,
                                           temp_count,dens_count:dens_count+len(gcr_vals)] = gcr_vals
            dens_count = dens_count + len(gcr_vals)
            if(dens_count == len(adf11['input_file']['dens'])):
                temp_count = temp_count + 1
                dens_count = 0
        gcr_line = f.readline()
