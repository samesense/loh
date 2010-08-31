"""Parse exome data to find minor allele frequency.
   Only look at calls w/ 8x coverage or more.
   Ignore strange chromosomes (including M).
   Output is in BED format 
   (http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED).
"""
from collections import defaultdict
import os, math

def get_freqs(call, str_length, str):
    """Count minor allele freq by matching . or ,
       Return difference between major and minor"""

    count = 0
    for s in str:
        if s == '.' or s == ',':
            count += 1
    # I think Murim's just taking the abs difference in freq counts
    return abs(float(2*count-str_length)/float(str_length))

BED_locations = {}
data_dir = 'data/exome/'
for afile in os.listdir(data_dir):
    if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
        file = os.path.join(data_dir, afile)
        sample2alleles = defaultdict(dict)
        with open(file) as f:
            for line in f:
                sp = line.split('\t')
#                print str(sp)
                samples = sp[8].split('-')
                idx = 17
                for s in samples:
                    chr = sp[2].split('chr')[1]
                    pos = sp[3]
                    quality = float(sp[idx+1])
                    coverage = int(sp[idx+2])
                    # use quality > 100
                    # keep if using at least 8x coverage
                    # and if the chr is standard
                    # i.e. no 6_cox_hap2
                    if coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
                        if 'XY' == chr:
                            chr = '100'
                        elif 'X' == chr:
                            chr = '98'
                        elif 'Y' == chr:
                            chr = '99'
                        sample2alleles[s][chr + '.' + pos] = get_freqs(sp[idx],
                                                                       coverage,
                                                                       sp[idx+3])
                        # keep track of genome locations
                        # for later conversion
                        BED_chr = sp[2]
                        BED_locations[BED_chr + ':' + pos + ':' + str(int(pos)+1)] = True
                    idx += 5
        outdir = os.path.join('working',
                              afile.replace('.', '_'))
        os.system('mkdir -p ' + outdir)
        for s in sample2alleles:
            ofile = os.path.join(outdir, 'hg19_chr_pos.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tMinFreq\n')
                for chr_pos in sample2alleles[s]:
                    chr, pos = chr_pos.split('.')
                    allele_freq = sample2alleles[s][chr_pos]
                    o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')
with open('working/BED_pos_hg19', 'w') as f:
    for position in BED_locations:
        chr, st, end = position.split(':')
        f.write(chr + '\t' + st + '\t' + end + '\n')
    
