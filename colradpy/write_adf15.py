import numpy as np
import itertools
from datetime import date

# group subset of an array
# this is so there are 8 entries per line
# in the adf15 file
# see stackexchangeG
def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = tuple(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk


def write_adf15(fil_name,pec_inds, wave, pecs, pec_lvls, element_sym,
                charge_state, dens, temp, metas, ion_pot, user='', atomic='', num=8,
                wave_lims=None):
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

    pecs[pecs<=1e-40] = 1e-40 # sets effective zero

    # Filter PECs by strength
    tol = 1e-17
    pec_inds_fil = np.array([], dtype='int')
    for pp in pec_inds:
        if np.mean(pecs[pp, :,:,:]) >= tol:
            pec_inds_fil = np.append(pec_inds_fil, pp)

    # Further filters PECs by wavelength
    if wave_lims is not None:
        pec_inds_fil2 = np.where(
            (wave >= wave_lims[0])
            & (wave <= wave_lims[1])
            )[0]

        # Takes intersection
        pec_inds_fil = list(set(pec_inds_fil) & set(pec_inds_fil2))

    pecs = pecs[pec_inds_fil]#only take PECs requested, dim(ntrans,nmeta,ntemp,ndens)
    wave = wave[pec_inds_fil]#only take wavelengths requested
    pec_lvls = pec_lvls[pec_inds_fil,:]#only take PEC levels requested

    #first line in the adf15 file format
    # nothing really important here as its all redunant information later in file
    first_line = '   ' + str(len(pec_inds_fil)*pecs.shape[1]) + '    /' + element_sym + \
                 ' ' + str(charge_state) + ' PHOTON EMISSIVITY COEFFICIENTS/\n'

    f.write(first_line)

    # Array to store transition info
    trans = []
    isel = 0

    for j in range(0,np.shape(pecs)[1]):#loop over metastables
        for i in range(0, len(wave)):#loop over PECs
            #letting the user know if the PEC is excitation or recombination
            if(j <len(metas)):
                typ = 'EXCIT' 
            elif(j < len(metas) + len(ion_pot)):
                 typ = 'RECOM'
            else:
                 typ = 'CHEXC'

            # Increases isel index
            isel += 1

            # Transition infor
            upr = pec_lvls[i,0]
            lwr = pec_lvls[i,1]
            trans.append(
                'C'+str(isel).rjust(6, ' ')
                + "{:1.5F}".format(wave[i]).rjust(15, ' ')
                + (
                    str(upr) 
                    +'('+str(atomic['S'][upr])+')'
                    +str(atomic['L'][upr])
                    +'('+str(atomic['w'][upr])+')'
                    ).rjust(17, ' ')
                +'-'
                + (
                    str(lwr) 
                    +'('+str(atomic['S'][lwr])+')'
                    +str(atomic['L'][lwr])
                    +'('+str(atomic['w'][lwr])+')'
                    ).rjust(15, ' ')
                + typ.rjust(8, ' ')
                + '\n'
                )

            #header information for the PEC, important things here are the wavelength
            # number of temps and dens and metastable level
            pec_header = '   ' + "{0:.5f}".format(wave[i]) + ' A  ' + str(len(dens)) +\
            '  ' + str(len(temp)) + ' /FILMEM = bottom  /TYPE = ' +typ +  ' /INDM = ' +\
            'T'+ ' /ISEL = ' + str(isel)+'\n'
            #str(j) + ' /ISEL = ' + str(isel)+'\n'

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
    f.write('C  NUCLEAR CHARGE =  '+str(element_sym)+'\n')
    f.write('C  ION CHARGE +1  =  '+str(charge_state+1)+'\n')
    f.write('C\n')
    f.write('C  SPECIFIC ION FILE  : '+user['file_loc']+'\n')
    f.write('C\n')
    f.write('C  TABULATION       : photon emissivity coefft (te,ne)\n')
    f.write('C  UNITS            : ph. emis coef[cm^3 s^-1]; te [eV]; ne [cm^-3]\n')
    f.write('C\n')    
    if(type(user) == dict):
        for i in user.keys():
            f.write('C '+ i + ' ' +str(user[i]) + "\n")

    if(type(atomic) == dict):
        f.write('C\n')
        f.write('C' + 'lvl'.rjust(5, ' ') +'config'.rjust(20,' ') +
                    'S'.rjust(5,' ')+
                    'L'.rjust(5,' ')+
                    'w'.rjust(5,' ')+
                    'energy [cm^-1]\n'.rjust(18, ' '))

        f.write('C-----------------------------------------------------------------------\n')
        f.write('C\n')        
        for ii in range(0,len(atomic['config'])):
            
            f.write('C' + str(ii).rjust(5, ' ') + str(atomic['config'][ii]).rjust(20, ' ')+\
                    str(atomic['S'][ii]).rjust(5, ' ') + \
                    str(atomic['L'][ii]).rjust(5, ' ') + \
                    str(atomic['w'][ii]).rjust(5, ' ') + \
                    "{:0.2F}".format(atomic['energy'][ii]).rjust(18, ' ')+'\n')

    # Transitions header
    f.write('C\n')
    f.write('C-----------------------------------------------------------------------\n')
    f.write('C\n')
    f.write(
        "C   ISEL     WVLEN [A]     TRANSITION lvl(S)L(w)           TYPE\n"
        )
    f.write(
        "C  -----  ------------  --------------------------------  -----\n"
        )
    for tt in np.arange(len(trans)):
        f.write(trans[tt])

    # Disclaimer
    today = date.today().strftime("%b-%d-%Y")
    f.write('C\n')
    f.write('C-----------------------------------------------------------------------\n')
    f.write('C\n')
    f.write('C\tPRODUCER:\t Conor Perks, cjperks@psfc.mit.edu\n')
    f.write('C\tDATE:    \t'+ today)
    f.write('\n')
    f.write('C\n')
    f.write('C-----------------------------------------------------------------------\n')
            
    f.close()
