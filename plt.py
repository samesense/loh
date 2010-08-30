from collections import defaultdict
import re

def get_freqs(call, str_length, str):
    """Count minor allele freq by matching . or ,"""

    count = 0
    for s in str:
        if s == '.' or s == ',':
            count += 1
    freq = min([float(count)/float(str_length),
                float(float(str_length)-count)/float(str_length)])
    print call, str_length, str, freq
    return freq

sample2alleles = defaultdict(dict)
with open('data/perry.txt') as f:
    for line in f:
        sp = line.split('\t')
        samples = sp[8].split('-')
        idx = 17
        for s in samples:
            chr = sp[2].split('chr')[1]
            pos = sp[3]
            # keep if using at least 8x coverage
            if int(sp[idx+2]) >= 8:
                sample2alleles[s][chr + '.' + pos] = get_freqs(sp[idx],
                                                               sp[idx+2],
                                                               sp[idx+3])
            idx += 5

for s in sample2alleles:
    ofile = 'working/chr_pos.' + s
    with open(ofile, 'w') as o:
        o.write('CHR\tMapInfo\tMinFreq\n')
        for chr_pos in sample2alleles[s]:
            chr, pos = chr_pos.split('.')
            allele_freq = sample2alleles[s][chr_pos]
            o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')
