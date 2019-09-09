import pymysql
from fractions import Fraction
import re




####################################################################################################
def sort_energy_key(dic):
    if(dic['energy'] ==''):#some entries in nist do not have energy ''
        dic['energy'] = -1 #set them to -1
    return float(dic['energy'])


def sort_energy(dic):
    dic.sort(key=sort_energy_key)#sorting by energy value entries '' get -1 to be removed
    return [x for x in dic if not (-1 == float(x.get('energy')))] #remove levels with no energy


def remove_conf_parentage(dic):
    for i in range(0,len(dic)):
        dic[i]['conf'] = re.sub(r'\([^)]*\).', '', dic[i]['conf'])
    return dic


def remove_senority(dic):
    for i in range(0,len(dic)):
        if(dic[i]['term'] == ''):
            dic[i]['term'] = '-1'
        else:
            dic[i]['term'] = dic[i]['term'][0:2]

    return dic


def get_nist_levels(element,charge):
    #connect to MYsql database 'levels'
    connection_levels = pymysql.connect(host='localhost',
                                        user='root',
                                        password='spectra',
                                        db='levels',
                                        charset='utf8mb4',
                                        cursorclass=pymysql.cursors.DictCursor)

    cur_levels = connection_levels.cursor()#lets you jump through the database

    cur_levels.execute(     #Make actual querery database based on parameters
        'SELECT * FROM ASD_Levels WHERE element= %s AND spectr_charge=%s', #(%s <- param_1, %s <- param_2)
            (element,charge))    # (param_1, param_2)
    
    return cur_levels.fetchall()



def get_nist_clean(element,charge):

    return remove_senority(remove_conf_parentage(sort_energy(get_nist_levels(element,charge))))
    
