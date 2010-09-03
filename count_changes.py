"""Count mutatios in somatic_mutations file.
   Used while deciding what a somatic change is.
"""
from collections import defaultdict

def init_zero():
    return 0

counts = defaultdict(init_zero)
seen = {}
with open('working/somatic_mutations') as f:
    for line in f:
        exome_type, normal, cancer, n_call, c_call = line.strip().split('\t')
        counts[n_call + ':' + c_call] += 1
for s in counts:
    n_call, c_call = s.split(':')
    print n_call + '\t' + c_call + '\t' + str(counts[s])
