"""Take novel results from get_novel_snps and summarize
   by patient"""
import sys
from collections import defaultdict

def init_sample():
    return {'missing_dbsnp':0,
            'missing_lifton':0,
            'missing_lifton_dbsnp':0,
            'not_novel':0,
            'total':0,
            'either':0}

def get_novel(afile):
    """Grab novel SNPs"""

    novel = {}
    with open(afile) as f:
        for line in f:
            sp = line.strip().split('\t')
            novel[sp[0] + ':' + sp[1]] = True
    return novel

def printit(data, sample):
    print('%s\t%d\t%d\t%d\t%d' %
          (sample, data[sample]['total'], data[sample]['missing_dbsnp'], 
           data[sample]['missing_lifton'], 
           data[sample]['missing_lifton_dbsnp']))
    
mu2a_input_file = sys.argv[1]
novel_dbsnp_file = sys.argv[2]
novel_lifton_file = sys.argv[3]

novel_dbsnp = get_novel(novel_dbsnp_file)
novel_lifton = get_novel(novel_lifton_file)

sample2counts = defaultdict(init_sample)
with open(mu2a_input_file) as f:
    for line in f:
        (sample, chr, pos, ref, mut) = line.strip().split('\t')
        chrpos = chr + ':' + pos
        sample2counts[sample]['total'] += 1
        if chrpos in novel_dbsnp and chrpos in novel_lifton:
             sample2counts[sample]['missing_lifton_dbsnp'] += 1
        if chrpos in novel_lifton:
            sample2counts[sample]['missing_lifton'] += 1
        if chrpos in novel_dbsnp:
            sample2counts[sample]['missing_dbsnp'] += 1
        # if chrpos in novel_dbsnp or chrpos in novel_lifton:
        #     sample2counts[sample]['either'] += 1
        # if chrpos in novel_lifton and chrpos not in novel_dbsnp:
        #     sample2counts[sample]['missing_lifton_only'] += 1
        # elif chrpos in novel_dbsnp and chrpos not in novel_lifton:
        #     sample2counts[sample]['missing_dbsnp_only'] += 1
        # elif chrpos in novel_dbsnp and chrpos in novel_lifton:
        #     sample2counts[sample]['missing_lifton_dbsnp'] += 1
        # else:
        #     sample2counts[sample]['not_novel'] += 1
print 'Patient\tTotal SNV\tdbSNP Novel\tLifton Novel\tBoth Novel'
for sample in sample2counts:
    printit(sample2counts, sample)

