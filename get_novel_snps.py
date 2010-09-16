"""Compare somatic mutation input with MU2A SNP calls
   to return a list of novel SNPs (SNPs not found by
   MU2A"""
import sys, mu2a

mutation_file = sys.argv[1]
mu2a_file = sys.argv[2]

mutations = {}
with open(mutation_file) as f:
    for line in f:
        (chr, pos, normal, cancer) = line.strip().split('\t')
        mutations[chr+':'+pos] = (normal, cancer)

mu2a_output = mu2a.MU2A(mu2a_file)
mu2a_snps = mu2a_output.get_snp_positions()

for chrpos in mutations:
    if chrpos not in mu2a_snps:
        chr, pos = chrpos.split(':')
        normal, cancer = mutations[chrpos]
        print('%s\t%s\t%s\t%s' %
              (chr, pos, 
               normal, cancer))
