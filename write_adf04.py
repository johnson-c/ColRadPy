from decimal import Decimal
#from read_adf04 import *
#from colradpy import colradpy

#a = colradpy('/home/curtis/spirit/wlike_mons11#w0.dat',[0])

f = open('adf04_out.txt','w')


def write_first_line():
    top_line = a['element'] + ' ' + str(a['charge_state']) + '        ' + a['iz0'] + '         ' + a['iz1'] + '       ' + str(a['ion_pot']) + '(' + a['ion_term'] + ')'
    f.write(top_line)
    f.close()



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
