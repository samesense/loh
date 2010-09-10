"""Yong's all_non_ref files have some problems calling AA:AB...
   This file locates those errors."""

def bad_call(mutation_type, norm_call, cancer_call):
    """Return true if the mutation type doesn't jive w/ the calls"""

    if mutation_type == 'AB:AA':
        pass
    elif mutation_type == 'AB:BB':
        pass
    elif mutation_type == 'BB:AB':
        pass
    return False

with open('working/total.new_mutants') as f:
    for line in f:
        chr, pos, mutation_type, norm_call, cancer_call = line.strip().split('\t')
        if bad_call(mutation_type, norm_call, cancer_call):
            print line.strip()
