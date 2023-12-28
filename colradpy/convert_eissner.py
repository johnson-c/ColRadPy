import re

def convert_eissner(eissner_config):
    """ converts an eissner notation string for the configuration to standard LS

    Args:
       :param eissner_config: string in eissner notation format
       :type eissner_config: string

    :returns: string of configuration in standard LS format
    """

    nl_pairs = {'1': '1s', '2': '2s', '3': '2p', '4': '3s', '5': '3p', '6': '3d', '7': '4s', '8': '4p', '9': '4d',
                'A': '4f', 'B': '5s', 'C': '5p', 'D': '5d', 'E': '5f', 'F': '5g', 'G': '6s', 'H': '6p', 'I': '6d',
                'J': '6f', 'K': '6g', 'L': '6h', 'M': '7s', 'N': '7p', 'O': '7d', 'P': '7f', 'Q': '7g', 'R': '7h',
                'S': '7i', 'T': '8s', 'U': '8p', 'V': '8d', 'W': '8f', 'X': '8g', 'Y': '8h', 'Z': '8i', 'a': '8k',
                'b': '9s', 'c': '9p', 'd': '9d', 'e': '9f', 'f': '9g', 'g': '9h', 'h': '9i', 'i': '9k', 'j': '9l',
                'k': '10s', 'l': '10p', 'm': '10d', 'n': '10f', 'o': '10g', 'p': '10h', 'q': '10i', 'r': '10k',
                's': '10l', 't': '10m', 'u': '11s', 'v': '11p', 'w': '11d', 'x': '11f', 'y': '11g', 'z': '11h'}



    elec  = {'1': '1', '2': '2', '3': '3', '4': '4', '5': '5', '6': '6', '7': '7', '8': '8', '9': '9',
                'A': '10', 'B': '11', 'C': '12', 'D': '13', 'E': '14', 'F': '15', 'G': '16', 'H': '17', 'I': '18'}


    if('0' == eissner_config[1]):#somethings ADAS goes deep into the closed shells,
                              #not sure whats going on with that just remove
        eissner_config = eissner_config[re.search('5',eissner_config).start():]


    eissner_config = re.sub('0','A',eissner_config)
        
    eissner_config = re.sub( '^5', '', eissner_config )#sometimes there is 5 at start

    ls_notat = ''


    for i in range(0,len(eissner_config),3):

        num_elec = int(elec[eissner_config[i]]) #electron number in the first space
        orb_num = nl_pairs[eissner_config[i+1]] #orbital is in the second space
                                               #Note 5 is a delimiter in the 3rd spot

        #since the configuration is already being modified from adf04 file
        #just go ahead and convert to be consistent with NIST

        '''
        if(num_elec >1):#remove single electron from configuration string consistent with NIST
            ls_notat = ls_notat + orb_num + str(num_elec)
        else:
        '''
        ls_notat = ls_notat + orb_num + str(num_elec)

        if(i <len(eissner_config)-2):#add in the . separation, consistent with NSIT
            ls_notat = ls_notat + '.'
    return ls_notat.upper()
