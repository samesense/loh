"""Compare SNP-chip and exome_seq cals"""
import sys, call_class, yusik_cmp_snpchip_exomeseq, utils, os, random
from collections import defaultdict

def get_snpchip_calls(afile):
    """Load rs# and allele calls. 
       These were contructed w/ snp.py"""

    rs2alleles = {}
    with open(afile) as f:
        for line in f:
            rs, gentrain, alleles = line.strip().split('\t')
            rs2alleles[rs] = (''.join(alleles.split('/')), float(100)*float(gentrain))
    return rs2alleles

def eval_exome(sample_name, snpchip_rs2alleles_normal, snpchip_rs2alleles_cancer):
    """Compare exome seq call to snp chip call"""

    bin = 10
    data_dir = os.path.join('data/all_non_ref_hg19/',
                            sample_name)
    mutation_calls = call_class.calls(data_dir, sample_name).data
    evals = {}
    for exome_type in mutation_calls:
        evals[exome_type] = {'consensus_quality':{}}
#                             'gentrain':{}}
                             #'min_both':{}}
        for chrpos in mutation_calls[exome_type]:
            if 'snp' in mutation_calls[exome_type][chrpos]:
                snp = mutation_calls[exome_type][chrpos]['snp']
                if snp in snpchip_rs2alleles_normal or snp in snpchip_rs2alleles_cancer:
                    if snp in snpchip_rs2alleles_normal:
                        chip_calls_normal, gentrain_normal = snpchip_rs2alleles_normal[snp]
                    else:
                        chip_calls_normal, gentrain_normal = ('NA', float(0))
                    if snp in snpchip_rs2alleles_cancer:
                        chip_calls_cancer, gentrain_cancer = snpchip_rs2alleles_cancer[snp]
                    else:
                        chip_calls_cancer, gentrain_cancer = ('NA', float(0))
                    normal_call = mutation_calls[exome_type][chrpos]['N']['call']
                    cancer_call = mutation_calls[exome_type][chrpos]['T']['call']
                    normal_quality = utils.my_round(mutation_calls[exome_type][chrpos]['N']['consensus_quality'], bin)
                    cancer_quality = utils.my_round(mutation_calls[exome_type][chrpos]['T']['consensus_quality'], bin)
                    # gentrain_quality_normal = utils.my_round(gentrain_normal, bin)
                    # gentrain_quality_cancer = utils.my_round(gentrain_cancer, bin)
                    # min_both_normal_quality = utils.my_round(min([gentrain, normal_quality]), bin)
                    # min_both_cancer_quality = utils.my_round(min([gentrain, cancer_quality]), bin)
                    for (a_quality, a_sample, 
                         a_quality_type, a_gentrain) in ((normal_quality, 'N', 'consensus_quality', gentrain_normal),
                                                         (cancer_quality, 'T', 'consensus_quality', gentrain_cancer)):
                                                                      #(gentrain_quality, 'N', 'gentrain'),
                                                                      #(gentrain_quality, 'T', 'gentrain')):
    #                                                                  (min_both_normal_quality, 'N', 'min_both'),
     #                                                                 (min_both_cancer_quality, 'T', 'min_both')):
                            if a_gentrain > float(90):
                                if a_sample not in evals[exome_type][a_quality_type]:
                                    evals[exome_type][a_quality_type][a_sample] = {}
                                if a_quality not in evals[exome_type][a_quality_type][a_sample]:
                                    evals[exome_type][a_quality_type][a_sample][a_quality] = [0,0]

                    # if 'N' not in evals[exome_type]['consensus_quality']:
                    #     evals[exome_type]['consensus_quality']['N'] = {}
                    # if normal_quality not in evals[exome_type]['consensus_quality']['N']:
                    #     evals[exome_type]['consensus_quality']['N'][normal_quality] = [0,0]

                    # if 'T' not in evals[exome_type]['consensus_quality']:
                    #     evals[exome_type]['consensus_quality']['T'] = {}
                    # if cancer_quality not in evals[exome_type]['T']:
                    #     evals[exome_type]['consensus_quality']['T'][cancer_quality] = [0,0]
                    res = yusik_cmp_snpchip_exomeseq.check_nuc(chip_calls_normal, 
                                                               normal_call, snp)
                    if res != 'NA' and gentrain_normal > float(90):
                        if res:
                            evals[exome_type]['consensus_quality']['N'][normal_quality][0] += 1
                        #evals[exome_type]['gentrain']['N'][gentrain_quality][0] += 1
                        #evals[exome_type]['min_both']['N'][min_both_normal_quality][0] += 1
                        else:
                            evals[exome_type]['consensus_quality']['N'][normal_quality][1] += 1
                        #evals[exome_type]['gentrain']['N'][gentrain_quality][1] += 1
                        #evals[exome_type]['min_both']['N'][min_both_normal_quality][1] += 1
                    res = yusik_cmp_snpchip_exomeseq.check_nuc(chip_calls_cancer, 
                                                               cancer_call, snp)
                    if res != 'NA' and gentrain_cancer > float(90):
                        if res:
                            evals[exome_type]['consensus_quality']['T'][cancer_quality][0] += 1
                        #evals[exome_type]['gentrain']['T'][gentrain_quality][0] += 1
                        #evals[exome_type]['min_both']['T'][min_both_cancer_quality][0] += 1
                        else:
                            evals[exome_type]['consensus_quality']['T'][cancer_quality][1] += 1
                        #evals[exome_type]['gentrain']['T'][gentrain_quality][1] += 1
                        #evals[exome_type]['min_both']['T'][min_both_cancer_quality][1] += 1

    tmpr = 'tmpr' + str(random.randint(0,1000))
    with open(tmpr, 'w') as f:
        f.write('Exome\tQuality_type\tSample\tQuality\tPercent_Match\n')
        for exome_type in evals:
            for quality_type in evals[exome_type]:
                for sample_type in evals[exome_type][quality_type]:
                    qualities = evals[exome_type][quality_type][sample_type].keys()
                    qualities.sort()
                    qualities.reverse()
                    sums = [0,0]
                    for q in qualities:
                        sums[0] += evals[exome_type][quality_type][sample_type][q][0]
                        sums[1] += evals[exome_type][quality_type][sample_type][q][1]
                        f.write('%s\t%s\t%s\t%.2f\t%.2f\n' %
                                (exome_type, quality_type, sample_type, q, 
                                 float(100)*float(sums[0])/float(sum(sums))))
    tmpR = 'tmprinput' + str(random.randint(0,1000))
    with open(tmpR, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write("png('plots/yuiri_exome_chip_cmp.png')\n")
        f.write("snpchip_exome_cmp('" + tmpr + "')\n")
        f.write('dev.off()\n')
        f.write('q()\n')
    os.system('R CMD BATCH --vanilla ' + tmpR + ' tmpLog')
    os.system('rm tmpLog ' + tmpR + ' ' + tmpr)
    
def main():
    """Entry"""

#    snpchip_file = sys.argv[1]
#    exomeseq_name = sys.argv[2]
    snp_chip_file_normal = 'working/yuiri_normal.call'
    snp_chip_file_cancer = 'working/yuiri_cancer.call'
    exomeseq_name = 'yuiri'
    snpchip_normal = get_snpchip_calls(snp_chip_file_normal)
    snpchip_cancer = get_snpchip_calls(snp_chip_file_cancer)
    eval_exome(exomeseq_name, snpchip_normal, snpchip_cancer)

if __name__ == '__main__':
    main()
    
