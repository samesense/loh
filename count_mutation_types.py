"""Count AA:AA, etc for new_mutants file made w/ mutations.py"""
from collections import defaultdict


def init_zero(): return 0

counts = defaultdict(init_zero)
with open('working/total.new_mutants') as f:
    for line in f:
        chr, pos, mutation_type, norm_call, cancer_call = line.strip().split('\t')
        counts[mutation_type] += 1
for mut in counts:
    print mut, counts[mut]
        
