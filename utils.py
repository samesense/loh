"""Functions called multiple times"""

def init_zero(): return 0

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
