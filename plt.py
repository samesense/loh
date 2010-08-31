"""Parse exome data to find minor allele frequency.
   Only look at calls w/ 8x coverage or more.
   Ignore strange chromosomes (including M).
"""
from collections import defaultdict
import os

def get_freqs(call, str_length, str):
    """Count minor allele freq by matching . or ,"""

    count = 0
    for s in str:
        if s == '.' or s == ',':
            count += 1
    freq = min([float(count)/float(str_length),
                float(float(str_length)-count)/float(str_length)])
    #print call, str_length, str, freq
    return freq

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
                    # keep if using at least 8x coverage
                    # and if the chr is standard
                    # i.e. no 6_cox_hap2
                    if int(sp[idx+2]) >= 8 and '_' not in chr and 'M' not in chr:
                        if 'XY' == chr:
                            chr = '100'
                        elif 'X' == chr:
                            chr = '98'
                        elif 'Y' == chr:
                            chr = '99'
                        sample2alleles[s][chr + '.' + pos] = get_freqs(sp[idx],
                                                                       sp[idx+2],
                                                                       sp[idx+3])
                    idx += 5
        outdir = os.path.join('working',
                              afile.replace('.', '_'))
        os.system('mkdir -p ' + outdir)
        for s in sample2alleles:
            ofile = os.path.join(outdir, 'chr_pos.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tMinFreq\n')
                for chr_pos in sample2alleles[s]:
                    chr, pos = chr_pos.split('.')
                    allele_freq = sample2alleles[s][chr_pos]
                    o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')
