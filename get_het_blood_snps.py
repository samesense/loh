"""Extract heterozygous SNPs from blood data.
   Require that a SNP is heterozygous in all
   normal samples."""
import os
from collections import defaultdict

normals = ('yuiriN', 'yuakerN', 'yuiskiaN', 'yusanN',
           'yunoca091283P')

hets = defaultdict(dict)
data_dir = 'data/exome/'
for afile in os.listdir(data_dir):
    if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
        file = os.path.join(data_dir, afile)
        sample2alleles = defaultdict(dict)
        with open(file) as f:
            for line in f:
                sp = line.split('\t')
                samples = sp[8].split('-')
                idx = 17
                for s in samples:
                    if s in normals:
                        chr = sp[2].split('chr')[1]
                        pos = sp[3]
                        quality = float(sp[idx+1])
                        coverage = int(sp[idx+2])
                        call = sp[idx]
                        # use quality > 100
                        # keep if using at least 8x coverage
                        # and if the chr is standard
                        # i.e. no 6_cox_hap2
                        if call not in ('A', 'T', 'G', 'C') and coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
                            if 'XY' == chr:
                                chr = '100'
                            elif 'X' == chr:
                                chr = '98'
                            elif 'Y' == chr:
                                chr = '99'
                            hets[chr+':'+pos][s] = True
                    idx += 5

# print out SNPs that are heterozygous in all samples
for chr_pos in hets:
    #if len(hets[chr_pos]) == len(normals):
    chr, pos = chr_pos.split(':')
    print chr + '\t' + pos
