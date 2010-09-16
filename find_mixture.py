"""This is the cumulation of the CNV/LOH work. We need the mixture of normal/cancer cells to help call somatic mutations. I will find the mixture per chromosome by looking at copy loss regions, as indicated by CNV-seq, and then look at the allele differences (Murim's plot) to determine the mixture by looking at the difference between 0.5 and the observed normal/cancer allele frequency differences."""
import sys
from collections import defaultdict

def avg_diffs_per_chr(allele_diffs, cnvs):
    """For each chr, use the lost of copy CNVs to compuete the averag allele diff per chr"""

    for chr in allele_diffs:
        sum_diffs = float(0)
        total_calls = 0
        c = 'chr' + chr
        if c in cnvs:
            for pos in allele_diffs[chr]:
                for cnv_st, cnv_end in cnvs[c]:
                    if pos < cnv_st:
                        break
                    if pos >= cnv_st and pos <= cnv_end:
                        sum_diffs += allele_diffs[chr][pos]
                        total_calls += 1
                        break
        if total_calls:
            avg_diff = str(sum_diffs/float(total_calls))
            print('%s\t%d\t%s' %
                  (chr, total_calls, avg_diff))
            
def get_cnvs(afile):
    """Load los of copy CNVs computed by cnv-seq"""

    cnvs = defaultdict(list)
    with open(afile) as f:
        f.readline()
        for line in f:
            sp = line.strip().split('\t')
            if len(sp) == 7:
                cnv, chr, st, end, size, log2, pval = line.strip().split('\t')
                if log2 != 'NA':
                    # check for copy loss
                    if float(log2) < float(-.8):
                        cnvs[chr].append((int(st), int(end)))
    for chr in cnvs:
        cnvs[chr].sort()
    return cnvs

def get_allele_diffs(afile):
    """Load Murim's normal/cancer diffs for heterozygous calls in normal."""

    allele_diffs = defaultdict(dict)
    with open(afile) as f:
        f.readline()
        for line in f:
            chr, pos, diff = line.strip().split('\t')
            allele_diffs[chr][int(pos)] = float(diff)
    return allele_diffs

def main():
    """Load CNVs and Murim allele diffs and find mixture per chromosome."""

    murim_diff_file = sys.argv[1] # working/exome_aa_chg/hg19_murim.yusan
    cnv_file = sys.argv[2] # working/cnv_seq/CNV/exome.aa_chg.yusan

    allele_diffs = get_allele_diffs(murim_diff_file)
    cnvs = get_cnvs(cnv_file)
    avg_diffs_per_chr(allele_diffs, cnvs)

if __name__ == '__main__':
    main()
