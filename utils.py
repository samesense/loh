"""Functions called multiple times"""

def check_input(afile):
    """See if there is data in here"""

    line_counts = 0
    with open(afile) as f:
        for line in f:
            line_counts += 1
    if line_counts > 10:
        return True
    else:
        return False
