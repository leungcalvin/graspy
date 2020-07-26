import numpy as np
import pandas as pd
import re

# fname = 'rlevelseV.out'
# fname = 'rlevelseV2.out'

def read_rlevelsout(fname):

    ## load all lines
    with open(fname, 'r', encoding='utf-8') as f:
        fstr = f.read()
    fstr = fstr.split('\n')

    N_line = len(fstr)

    # TODO: check if the right file

    ## read
    # find start and end line for data
    str_header = ' No Pos  J'
    for li in range(N_line):
        if fstr[li].find(str_header) != -1:
            li_header = li
            break

    li_unit = li_header + 1
    li_start = li_header + 3

    str_end = '----'
    for li in range(li_start, N_line):
        if fstr[li].find(str_end) != -1:
            li_end = li - 1
            break

    # read energy unit (Levels, Splitting)
    mo = re.search('(cm\^-1|eV)',fstr[li_unit])
    unit_E = mo.group()

    # read out data
    pattern = '^\s*\d+\s+\d+' \
              '\s+(?P<J>\d+(/\d+)?)' \
              '\s+(?P<Parity>(\+|-)?)' \
              '\s+(?P<E_tot>-?\d+\.\d+)' \
              '(\s+(?P<Level>-?\d+\.\d+))?' \
              '(\s+(?P<Splitting>-?\d+\.\d+))?' \
              '(\s+(?P<Configuration>[^\s].+[^\s]))?' \
              '\s*$'
    # pattern = '^\s*\d+\s+\d+' \
    #           '\s+(?P<J>\d+(/\d+)?)' \
    #           '\s+(?P<Parity>(\+|-)?)'
    reo = re.compile(pattern)

    LD = [] # list of dicts
    for li in range(li_start, li_end + 1):
        mo = reo.search(fstr[li])
        LD.append(mo.groupdict(default=np.nan))

    N_level = li_end - li_start + 1

    # convert string to number
    label_2num = ['E_tot', 'Level', 'Splitting']
    for i in range(N_level):
        for label in label_2num:
            LD[i][label] = float(LD[i][label])


    ## put data in pandas.DataFrame
    # add unit
    D_unit = {'J': '',
              'Parity': '',
              'E_tot': 'Hartree',
              'Level': unit_E,
              'Splitting': unit_E,
              'Configuration': ''}

    df = pd.DataFrame([D_unit] + LD, index=['Unit'] + list(range(N_level)))

    return df




