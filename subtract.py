import sys

def load(afile):
    d = {}
    with open(afile) as f:
        for line in f:
            d[line.strip()] = True
    return d
d1 = load(sys.argv[1])
d2 = load(sys.argv[2])
for k in d1:
    if k not in d2:
        print k
