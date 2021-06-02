import numpy as np
import itertools


# group subset of an array
# this is so there are 8 entries per line
# in the adf15 file
# see stackexchange
def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk


def write_adf15(fil_name,pec_inds, wave, pecs, element_sym,
                charge_state, dens, temp, metas, ion_pot, user='', atomic='', num=8):
    """Write all or a subset of PECs to an adf15 formatted file.
    standard adf15 units
    wave in [A]
    PEC in [ph cm3 s-1]
    dens in  [cm-3]
    temp in [eV]

    :param fil_name
    :type str

    :param pec_inds
    :type int arr

    :param wave
    :type float arr

    :param pecs
    :type float arr

    :param element_sym
    :type str

    :param charge_state
    :type int

    :param dens
    :type float arr

    :param temp
    :type float arr

    :param metas
    :type float_arr

    :param ion_pot
    :type float arr


    :param user
    :type dict

    :param atomic
    :type dict

    :param num
    :type int

    """

    f = open(fil_name,'w+')#open the file

    pecs = pecs[pec_inds]#only take PECs requested
    wave = wave[pec_inds]#only take wavelengths requested

    #first line in the adf15 file format
    # nothing really important here as its all redunant information later in file
    first_line = '   ' + str(len(wave)) + '    /' + element_sym + \
                 ' ' + str(charge_state) + ' PHOTON EMISSIVITY COEFFICIENTS/\n'

    f.write(first_line)

    for j in range(0,np.shape(pecs)[1]):#loop over metastables
        for i in range(0, len(wave)):#loop over PECs
            #letting the user know if the PEC is excitation or recombination
            if(j <len(metas)):
                typ = 'EXCIT' 
            elif(j < len(metas) + len(ion_pot)):
                 typ = 'RECOM'
            else:
                 typ = 'CHEXC'

            #header information for the PEC, important things here are the wavelength
            # number of temps and dens and metastable level
            pec_header = '   ' + "{0:.1f}".format(wave[i]*10) + ' A  ' + str(len(dens)) +\
            '  ' + str(len(temp)) + ' /FILMEM = bottom  /TYPE = ' +typ +  ' /INDM = ' +\
            str(j) + ' ISEL = ' + str(i+1)+'\n'

            f.write(pec_header)

            #write the density grid to the file 'num' entries per line
            for chunk in grouper(num,dens):
                f.write(' '.join(map('{:.2E}'.format,chunk)))
                f.write('\n')
            #write the temperature grid to the file 'num' entries per line            
            for chunk in grouper(num,temp):
                f.write(' '.join(map('{:.2E}'.format,chunk)))
                f.write('\n')
            #write the PEC grid to the file 'num' entries per line
            #all the densities then temperatures
            for d in range(0,len(dens)):
                for chunk in grouper(num,pecs[i,j,:,d]):
                    f.write(' '.join(map('{:.2E}'.format,chunk)))
                    f.write('\n')

    #write out calculation details
    f.write('C-----------------------------------------------------------------------\n')
    f.write('C\n')
    f.write('C  PHOTON EMISSIVITY COEFFICIENTS:\n')
    f.write('C\n')
    f.write('C\n')
    f.write('C  INFORMATION\n')
    f.write('C  -----------\n')
    f.write('C\n')
    f.write('C   ***MADE WITH ColRadPy***\n')
    f.write('C  NUCLEAR CHARGE =  4\n')
    f.write('C  ION CHARGE +1  =  3\n')
    f.write('C\n')
    f.write('C  SPECIFIC ION FILE  :\n')
    f.write('C\n')
    f.write('C\n')    
    if(type(user) == dict):
        for i in user.keys():
            f.write('C '+ i + ' ' +str(user[i]) + "\n")

    if(type(atomic) == dict):
        f.write('C  ' + 'lvl' +'  '+ 'config' +
                    '  ' + 'S'+
                    '  ' + 'L'+
                    '  ' + 'w'+
                    '  ' + 'eng\n')

        f.write('C-----------------------------------------------------------------------\n')
        f.write('C\n')        
        for ii in range(0,len(atomic['config'])):
            
            f.write('C  ' + str(ii) +'  '+ str(atomic['config'][ii])+\
                    '  ' + str(atomic['S'][ii]) + \
                    '  ' + str(atomic['L'][ii]) + \
                    '  ' + str(atomic['w'][ii]) + \
                    '  ' + str(atomic['energy'][ii])+'\n')
            
    f.close()
