import urllib.request
import numpy as np


url = "http://open.adas.ac.uk/download/adf04/"
h_arr = np.array(['[adas][1/ha00_ls][h0.dat','adas][1/ha00_ls][h1.dat'])
h_file = np.array(['h0_adf04','h1_adf04'])

he_arr = np.array(['adas][2/mom97_ls][he0.dat','adas][2/mom97_ls][he1.dat'])
he_file = np.array(['he0_adf04','he1_adf04'])

li_arr = np.array(['adas][3/cpb02_ls][li0.dat','adas][3/cpb02_ls][li1.dat','adas][3/cpb02_ls][li2.dat'])
li_file = np.array(['li0_adf04','li1_adf04','li2_adf04'])

be_arr = np.array(['adas][4/cpb03_ls][be0.dat','adas][4/cpb03_ls][be1.dat','adas][4/cpb03_ls][be2.dat','adas][4/cpb03_ls][be3.dat'])
be_file = np.array(['be0_adf04','be1_adf04','be2_adf04','be3_adf04'])
                  

be_arr = np.array(['adas][4/cpb03_ls][be0.dat','adas][4/cpb03_ls][be1.dat','adas][4/cpb03_ls][be2.dat','adas][4/cpb03_ls][be3.dat'])
be_file = np.array(['be0_adf04','be1_adf04','be2_adf04','be3_adf04'])

b_arr = np.array(['cophps][b/dw/ls][b0.dat','copaw][be/belike_lfm14][b1.dat','copaw][li/lilike_lgy10][b2.dat','cophps][he/dw/ls][b3.dat','hlike/hlike_cpb02][b4.dat'])
b_file = np.array(['b0_adf04','b1_adf04','b2_adf04','b3_adf04'])






urllib.request.urlretrieve(url+h_arr[1], 'h1_adf04') 




def fix_whitespace(fname):
    """ Fix whitespace in a file """
    with open(fname, "rb") as fo:
        original_contents = fo.read()
    # "rU" Universal line endings to Unix
    with open(fname, "rU") as fo:
        contents = fo.read()
    lines = contents.split("\n")
    fixed = 0
    for k, line in enumerate(lines):
        new_line = line.rstrip()
        if len(line) != len(new_line):
            lines[k] = new_line
            fixed += 1
    with open(fname, "wb") as fo:
        fo.write("\n".join(lines))
    if fixed or contents != original_contents:
        print("************* %s" % os.path.basename(fname))
    if fixed:
        slines = "lines" if fixed > 1 else "line"
        print("Fixed trailing whitespace on %d %s" \
              % (fixed, slines))
    if contents != original_contents:
        print("Fixed line endings to Unix (\\n)")
