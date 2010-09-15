"""Tools for dealing w/ bed tools"""
from collections import defaultdict
import bisect

def load_bed(afile, track):
    """Read BED file as {} of chr2st2end and {} of pos lists for each chr"""

    found_track = False
    with open(afile) as f:
        for line in f:
            if 'track' in line:
                if track in line:
                    chr2st2end = defaultdict(dict)
                    chr2posLs = defaultdict(list)
                    found_track = True
                elif found_track:
                    return (chr2st2end, chr2posLs)
            elif found_track:
                chr, st, end = line.strip().split('\t')
                start = int(st)
                chr2st2end[chr][start] = int(end)
                chr2posLs[chr].append(start)
    if found_track:
        return (chr2st2end, chr2posLs)
    else:
        sys.stderr.write(track + 'not in ' + afile)

def check_interval(pos_to_check, start, end):
    """See if pos_to_check is in this interval"""

    if pos_to_check >= start and pos_to_check <= end:
        return True
    return False

def find_ge(a, x):
    """Find leftmost item greater than or equal to x"""

    i = bisect.bisect_left(a, x)
    if i == len(a):
        i -= 1
    #if i != len(a):
    return (a[i], i)
    #raise ValueError

def find_location_in_bed(chr, pos, bed_struct, chr2st2end):
    """Find this position in bed_struct. If found, return st and stp, 
       else False"""

    if 'chr' not in chr:
        raise ValueError
    
    left_most_ge, idx = find_ge(bed_struct[chr],
                                pos)
    if left_most_ge == pos:
        start = bed_struct[chr][idx]
        return left_most_ge
    else:
        start = bed_struct[chr][idx-1]

    if check_interval(pos, start, chr2st2end[chr][start]):
        return start
    else:
        return False
