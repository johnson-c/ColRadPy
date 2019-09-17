import numpy as np

def get_nist_txt(el,ch):

    tmp = np.genfromtxt('../atomic/nist_energies/'+el.lower()+'_'+str(ch),dtype=None,delimiter=',')
    levels = []
    for i in range(0,len(tmp)):
        tmp_dic = {}
        tmp_dic['conf']   = tmp[i][0].decode('UTF-8')
        tmp_dic['term']   = tmp[i][1].decode('UTF-8')
        if(type(tmp[i][2]) == np.bytes_):
            tmp_dic['j_val']  = tmp[i][2].decode('UTF-8')
        else:
            tmp_dic['j_val'] = str(tmp[i][2])
        tmp_dic['energy'] = tmp[i][3]
        levels.append(tmp_dic)

    return levels
