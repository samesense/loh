"""Recreate Murim's plot for paired normal/cancer samples."""
import os, sys, global_settings
from collections import defaultdict

bases = ('A', 'T', 'G', 'C')

# cancer, normal
pairs = global_settings.pairs
cancer2normal = {}
for cancer, normal in pairs:
    cancer2normal[cancer] = normal
normals = [x[1] for x in pairs]
cancers = [x[0] for x in pairs]
all = []
for n in normals: all.append(n)
for c in cancers: all.append(c)

def init_zero():
    return 0

def get_normal_val(ref_allele, str_length, str):
    """Count minor allele frequency."""

    ref_count = 0
    base_counts = defaultdict(init_zero)
    for s in str:
        if s == '.' or s == ',':
            ref_count += 1
        elif s.upper() in bases:
            base_counts[s.upper()] += 1
    base_count_ls = []
    for base in base_counts:
        base_count_ls.append((base_counts[base], base))
    base_count_ls.sort()
    if ref_allele == base_count_ls[-1][1]:
        use_base = base_count_ls[-2][1]
    else:
        use_base = base_count_ls[-1][1]
        
    return (float(base_counts[use_base])/float(str_length), use_base)

def get_het_normal(file):
    """Grab hetero SNPs from normal samples"""
    
    normal2het = defaultdict(dict)
    with open(file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            for s in samples:
                if s in normals:
                    ref_allele = sp[10]
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
                        normal2het[s][chr+':'+pos] = get_normal_val(ref_allele,
                                                                    coverage,
                                                                    sp[idx+3])
                idx += 5
    return normal2het

def get_norm_sample(sample):
    """If this sample is cancer, and pairs
       with a normal, return that normal.
       Otherwise return False."""

    if sample in cancer2normal:
        return cancer2normal[cancer]
    else:
        return False

def get_freq(ref_allele, call, str_length, str, normal_freq_and_base):
    """Take diff btwn normal_base_freq and cancer"""

    ref_count = 0
    base_counts = defaultdict(init_zero)
    for s in str:
        if s == '.' or s == ',':
            ref_count += 1
        elif s.upper() in bases:
            base_counts[s.upper()] += 1
    normal_freq, normal_base = normal_freq_and_base
    if normal_base in base_counts:
        return abs(normal_freq - float(base_counts[normal_base])/float(str_length))
    elif normal_base == call:
        return abs(normal_freq - float(ref_count)/float(str_length))
    else:
        return float(0)

data_dir = 'data/exome/'
for afile in os.listdir(data_dir):
    if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
        file = os.path.join(data_dir, afile)
        normal2het = get_het_normal(file)
        sample2alleles = defaultdict(dict)
        with open(file) as f:
            for line in f:
                sp = line.split('\t')
                samples = sp[8].split('-')
                idx = 17
                for s in samples:
                    norm_sample = get_norm_sample(s)
                    if norm_sample:
                        chr = sp[2].split('chr')[1]
                        pos = sp[3]
                        quality = float(sp[idx+1])
                        coverage = int(sp[idx+2])
                        ref_allele = sp[10]
                        # use quality > 100
                        # keep if using at least 8x coverage
                        # and if the chr is standard
                        # i.e. no 6_cox_hap2
                        # sample must be het in normal or paired w/ het in normal
                        if chr+':'+pos in normal2het[norm_sample] and coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
                            if 'XY' == chr:
                                chr = '100'
                            elif 'X' == chr:
                                chr = '98'
                            elif 'Y' == chr:
                                chr = '99'
                            sample2alleles[s][chr + ':' + pos] = get_freq(ref_allele,
                                                                          sp[idx],
                                                                          coverage,
                                                                          sp[idx+3],
                                                                          normal2het[norm_sample][chr+':'+pos])
                    idx += 5
        outdir = os.path.join('working',
                              afile.replace('.', '_'))
        os.system('mkdir -p ' + outdir)
        for s in sample2alleles:
            ofile = os.path.join(outdir, 'hg19_murim.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tDiff\n')
                for chr_pos in sample2alleles[s]:
                    chr, pos = chr_pos.split(':')
                    allele_freq = sample2alleles[s][chr_pos]
                    o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')
        for s in normal2het:
            ofile = os.path.join(outdir, 'hg19_murim.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tDiff\n')
                for chr_pos in normal2het[s]:
                    chr, pos = chr_pos.split(':')
                    freq, base = normal2het[s][chr_pos]
                    allele_freq = abs(float(.5) - freq)
                    o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')




