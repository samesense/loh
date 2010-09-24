"""Compare SNP-chip and exome_seq cals"""
import sys, call_class, yusik_cmp_snpchip_exomeseq, utils, os
from collections import defaultdict

def get_snpchip_calls(afile):
    """Load rs# and allele calls. 
       These were contructed w/ snp.py"""

    rs2alleles = {}
    with open(afile) as f:
        for line in f:
            rs, alleles = line.strip().split('\t')
            rs2alleles[rs] = ''.join(alleles.split('/'))
    return rs2alleles

def eval_exome(sample_name, snpchip_rs2alleles):
    """Compare exome seq call to snp chip call"""

    data_dir = os.path.join('data/all_non_ref_hg19/',
                            sample_name)
    mutation_calls = call_class.calls(data_dir, sample_name).data
    evals = {}
    for exome_type in mutation_calls:
        evals[exome_type] = {}
        for chrpos in mutation_calls[exome_type]:
            if 'snp' in mutation_calls[exome_type][chrpos]:
                snp = mutation_calls[exome_type][chrpos]['snp']
                if snp in snpchip_rs2alleles:
                    chip_calls = snpchip_rs2alleles[snp]
                    normal_call = mutation_calls[exome_type][chrpos]['N']['call']
                    cancer_call = mutation_calls[exome_type][chrpos]['T']['call']
                    normal_quality = round(mutation_calls[exome_type][chrpos]['N']['consensus_quality'], -1)
                    cancer_quality = round(mutation_calls[exome_type][chrpos]['T']['consensus_quality'], -1)
                    if 'N' not in evals[exome_type]:
                        evals[exome_type]['N'] = {}
                    if normal_quality not in evals[exome_type]['N']:
                        evals[exome_type]['N'][normal_quality] = [0,0]

                    if 'T' not in evals[exome_type]:
                        evals[exome_type]['T'] = {}
                    if cancer_quality not in evals[exome_type]['T']:
                        evals[exome_type]['T'][cancer_quality] = [0,0]

                    if yusik_cmp_snpchip_exomeseq.check_nuc(chip_calls, normal_call, snp):
                        evals[exome_type]['N'][normal_quality][0] += 1
                    else:
                        evals[exome_type]['N'][normal_quality][1] += 1

                    if yusik_cmp_snpchip_exomeseq.check_nuc(chip_calls, cancer_call, snp):
                        evals[exome_type]['T'][cancer_quality][0] += 1
                    else:
                        evals[exome_type]['T'][cancer_quality][1] += 1
    for exome_type in evals:
        for sample_type in evals[exome_type]:
            qualities = evals[exome_type][sample_type].keys()
            qualities.sort()
            qualities.reverse()
            sums = [0,0]
            for q in qualities:
                sums[0] += evals[exome_type][sample_type][q][0]
                sums[1] += evals[exome_type][sample_type][q][1]
                print exome_type, sample_type, q, sums[0], sums[1], float(100)*float(sums[0])/float(sum(sums))

def main():
    """Entry"""

#    snpchip_file = sys.argv[1]
#    exomeseq_name = sys.argv[2]
    snp_chip_file = 'working/yuiri_normal.call'
    exomeseq_name = 'yuiri'
    snpchip = get_snpchip_calls(snp_chip_file)
    eval_exome(exomeseq_name, snpchip)

if __name__ == '__main__':
    main()
    
