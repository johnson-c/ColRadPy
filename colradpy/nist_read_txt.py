import re


#J of a level, such as 2, 3/2 or 1/2 or 3/2. Unresolved levels list several J
#values separated by commas (3P,0,1,2), which continue the J column.
J_VALUE = re.compile(r'\d+(/\d+)?( or \d+(/\d+)?)*')


def _split_level(line):
    """Configuration, term, J and energy of one 'conf,term,J,energy' line.

       Commas can also be inside the configuration (2,3.sp- or 5s<1/2,F=4>.5g)
       or between the J values of an unresolved level. The energy is the last
       field, extra J values before it belong to J and the configuration takes
       the remaining fields before the term.
    """
    fields = line.split(',')
    j_start = len(fields) - 2
    while(j_start > 2 and J_VALUE.fullmatch(fields[j_start - 1])):
        j_start = j_start - 1
    return (','.join(fields[:j_start - 1]), fields[j_start - 1],
            ','.join(fields[j_start:-1]), float(fields[-1]))


def get_nist_txt(path,el,ch):

    levels = []
    with open(path+el.lower()+'_'+str(ch)) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if(line.strip() == ''):
                continue
            tmp_dic = {}
            tmp_dic['conf'], tmp_dic['term'], tmp_dic['j_val'], tmp_dic['energy'] = _split_level(line)
            levels.append(tmp_dic)

    return levels
