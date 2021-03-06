"""Functions called multiple times"""
from collections import defaultdict

def my_round(x, base):
    """Round number x to any base"""

    return int(base * round(float(x)/base))    

def flip_bases(bases):
    """Take the complement of these two bases. Order does not matter"""

    new_bases = ''
    for base in bases:
        new_bases = new_bases + global_settings.comp[base]
    return new_bases

def check_bases(bases, base1, base2):
    if (bases[0] == base1 and bases[1] == base2) or (bases[0] == base2 and bases[1] == base1):
        return True
    else:
        return False

def bases2het(bases):
    """Convert 2 bases to single char code"""

    if bases[0] == bases[1]:
        return bases[0]
    if check_bases(bases, 'A', 'C'):
        return 'M'
    elif check_bases(bases, 'A', 'G'):
        return 'R'
    elif check_bases(bases, 'A', 'T'):
        return 'W'
    elif check_bases(bases, 'C', 'G'):
        return 'S'
    elif check_bases(bases, 'C', 'T'):
        return 'Y'
    elif check_bases(bases, 'G', 'T'):
        return 'K'

def init_zero(): return 0

def init_zero_lists(): return [0,0]

def init_zero_paired_list(): 
    return {'T':defaultdict(init_zero_lists),
            'N':defaultdict(init_zero_lists)}

def check_input(afile):
    """See if there is data in this file"""

    line_counts = 0
    with open(afile) as f:
        for line in f:
            line_counts += 1
    if line_counts > 10:
        return True
    else:
        return False
