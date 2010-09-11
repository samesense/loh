"""Count AA:AA, etc for new_mutants file made w/ mutations.py"""
from collections import defaultdict
import sys

def init_zero(): return 0

counts = defaultdict(init_zero)
# file 'working/total.new_mutants'
with open(sys.argv[1]) as f:
    for line in f:
        chr, pos, mutation_type, norm_call, cancer_call = line.strip().split('\t')
        counts[mutation_type] += 1
for mut in counts:
    print mut, counts[mut]
        
