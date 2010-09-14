import sys

s1 = {}
with open(sys.argv[1]) as f:
    for line in f:
        s1[line.strip()] = True
s2 = {}
with open(sys.argv[2]) as f:
    for line in f:
        s2[line.strip()] = True
for item in s1:
    if item in s2:
        print item
