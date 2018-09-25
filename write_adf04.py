import numpy as np

def write_first_line(fil,element,charge_state,iz0,iz1,ion_pot,ion_term,append=False):
    top_line = element + '+ ' + str(charge_state) + '        ' \
               + str(iz0) + '         ' + str(iz1) 
    for i in range(0,len(ion_pot)):
        top_line = top_line + '       '+ str(ion_pot[i]) + '(' +\
                   str(ion_term[i]) + ')\n'
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')
    f.write(top_line)
    f.close()
    return 

def write_level_info(fil, config, S, L, w, energy, zpla, zpla1,append=True):
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')
    if(len(config) == len(S) == len(L) == len(w) ==  \
       len(energy) == len(zpla) == len(zpla1)):
        for i in range(0,len(config)):
            line = ' ' *(5-len(str(i+1)))
            line = line + str(i+1)# level number starts at 1
            line = line + ' ' + config[i]
            spaces = ' ' * (8+len(max(config,key=len))-len(config[i]))
            line = line  + spaces + '(' + str(S[i]) + ')' + \
                   str(L[i]) + '( ' + str(w[i]) + ')'
            spaces_1 = ' ' * (15 - len(format(energy[i],'.4f')))
            line = line +spaces_1+ format(energy[i],'.4f')
            for j in range(0,len(zpla[i])):

                line = line+ '    {' + str(zpla1[i][j]) + '}' + \
                       format(zpla[i][j],'.2f')
            line = line + '\n'
            f.write(line)
        f.write('   -1\n')
    else:
        print('All entries must have same length')

    f.close()

def write_temperature(fil,temp,append=True):
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')

    line = ' 1.00    3          '
    for i in range(0,len(temp)):
        line  = line + ' ' + format(temp[i],'.2e').replace('e','')
    line = line + '\n'
    f.write(line)
    f.close()

def write_col_rates(fil,col_transitions,a_val,col_excit,inf_engy=np.array([])\
                    ,append=True,reorder=True):
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')
    for i in range(0,len(col_transitions)):
        spaces_1 = 6 - len(str(col_transitions[i][0]))
        spaces_2 = 6 - len(str(col_transitions[i][1]))
        line = ' '*spaces_1 + str(col_transitions[i][0]) + ' '*spaces_2 +\
               str(col_transitions[i][1])+' ' +format(a_val[i],'.2e').\
               replace('e','')+ ' ' + np.array2string(col_excit[i],\
                formatter={'float_kind':lambda x: "%.2e" % x}).\
                replace('e', '').replace('[','').replace(']','').\
                replace('\n','')

        if(inf_engy.any()):
            line = line + ' ' + format(inf_engy[i],'.2e').replace('e','')+'\n'
            
        else:
            line = line + '\n'
        
        f.write(line)
    f.close()


def write_other_rates(fil,idd,transitions,rates,append=True):
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')
    for i in range(0,len(transitions)):
        spaces_1 = 5 - len(str(transitions[i][0]))
        spaces_2 = 5 - len(str(transitions[i][1]))
        line = idd+' '*spaces_1 + str(transitions[i][0]) + \
               ' '*spaces_2 +'+'+ str(transitions[i][1]) + '         '
        line = line + np.array2string(rates[i],\
                formatter={'float_kind':lambda x: "%.2e" % x}).\
                replace('e', '').replace('[','').replace(']','').\
                replace('\n','') + '\n'
        f.write(line)
    f.close()

def finish_and_comment(append=True):
    if(append):
        f = open(fil,'a+')
    else:
        f = open(fil,'w+')    
    line = '   -1\n   -1   -1\n'




'''
def write_adf04_user_temp(fil,adf04):
    write_first_line(fil,adf04['element'],adf04['charge_state'],adf04['iz0'],adf04['iz1'],adf04['ion_pot'],adf04['ion_term'],append=False):
    write_level_info(fil,adf04[' config'],adf04[' S'],adf04[' L'],adf04[' w'],adf04[' energy'],adf04[' zpla'],adf04[' zpla1'],append=True):
    write_temperature(fil,adf04['user_temp_grid'],append=True):    
'''

    
'''    
def write_col_excit_line(i):
    num_spaces = 5
    line_string = ''
    upper = str(a['col_transitions'][i,0])
    lower = str(a['col_transitions'][i,1])
    for j in range(0,num_spaces - len(upper)):
        line_string = line_string + ' '
    line_string = line_string + upper
    for j in range(0,num_spaces - len(lower)):
        line_string = line_string + ' '
    line_string = line_string + lower
    line_string = line_string + ' ' 
    av_string = '%.2E' % Decimal(a['a_val'][i])

    if( '-' not in av_string):
        av_string = av_string.replace('E+','+')
    elif('-' in av_string):
        av_string = av_string.replace('E-','-')

    line_string = line_string + av_string
    line_string = line_string + ' '
    for j in range(0,len(a['col_excit'][i])):
        tran = '%.2E' % Decimal(a['col_excit'][i][j])

        if( '-' not in tran):
            tran = tran.replace('E+','+')
        elif('-' in tran):
            #tran = tran[1:]
            tran = tran.replace('E-','-')
        line_string = line_string + tran
        line_string = line_string + ' '
    line_string = line_string[0:len(line_string)-1]
    line_string = line_string +'\n'
    f.write(line_string)



for i in range(0,len(a['col_transitions'])):
    if(a['col_transitions'][i,0] < 441):
        write_col_excit_line(i)
f.close()
'''
