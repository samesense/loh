"""For Yusik, I have nucleotides for cancer and normal SNP chips,
   so I can compare the calls"""
import utils, global_settings
from collections import defaultdict

yusik_snpchip_file = 'data/snp_chip/YUSIK_SNPs_from_Chip.txt'
yusik_exomeseq_file = 'data/all_non_ref_hg19/yusik/exome.intron'

def check_nuc(snpchip, exome, snp):
    """Do snpchip and exome match?"""

    snp_base = utils.bases2het(snpchip)
    if snp_base == exome:
        return True
    else:
        if (snp_base == 'M' and exome == 'K') or (snp_base == 'K' and exome == 'M'):
            return True
        elif (snp_base == 'Y' and exome == 'R') or (snp_base == 'R' and exome == 'Y'):
            return True
        # elif (snp_base in ('A', 'C', 'T', 'G')) and (exome in ('A', 'C', 'T', 'G')):
        #     if global_settings.comp[snp_base] == exome:
        #         print 'flip', snp_base, exome, snp
        #         return True
        #     else:
        #         return False
        else:
            #print snp_base, exome
            return False

def cmp_snpchip(cancer_exome_snps, normal_exome_snps, flip):
    """Count right/wrong calls. Flip snpchip bases as needed.
       1 in flip means the bases need to be reversed."""

    cancer_eval = [0,0]
    normal_eval = [0,0]
    missing = {}
    with open(yusik_snpchip_file) as f:
        for line in f:
            sp = line.strip().split('\t')
            snp = sp[1]
            if snp in flip:
                cancer_call = sp[15]
                normal_call = sp[22]
                if snp in cancer_exome_snps:
                    if cancer_call != '--':
                        if flip[snp] == '1':
                            cancer_call = utils.flip_bases(cancer_call)
                        if check_nuc(cancer_call, cancer_exome_snps[snp], snp):
                            cancer_eval[0] += 1
                        else:
                            cancer_eval[1] += 1
                if snp in normal_exome_snps:
                    if normal_call != '--':
                        if flip[snp] == '1':
                            normal_call = utils.flip_bases(normal_call)
                        if check_nuc(normal_call, normal_exome_snps[snp], snp):
                            normal_eval[0] += 1
                        else:
                            normal_eval[1] += 1
            elif 'NC' not in line:
                missing[snp] = True
    print 'cancer/normal, hit/miss'
    print cancer_eval, float(cancer_eval[0]/float(sum(cancer_eval))), normal_eval, float(100)*float(normal_eval[0])/float(sum(normal_eval))
    print 'missing', len(missing)
    #print missing

def load_exome_snps(quality_cut, coverage_cut):
    """Get cancer/normal calls from Yong's file"""

    snps = defaultdict(dict)
    with open(yusik_exomeseq_file) as f:
        for line in f:
            sp = line.strip().split('\t')
            if int(sp[-2]) >= coverage_cut and float(sp[-8]) > quality_cut:
                snp = sp[5]
                if 'rs' in snp:
                    snps[sp[8]][snp] = sp[17]
    return snps

def check_dbsnp_illumina():
    """See if I need to flip the bases for these rs#s.
       Ensure each ss# maps to only one rs#"""

    rs2ss = {}
    problems = {}
    with open('data/dbsnp/ILLUMINA_Human_1M') as f:
        for line in f:
            if line[0:2] == 'ss':
                sp = line.strip().split('\t')
                ss = sp[0]
                rs = sp[4]
                flip_bases = sp[5]
                assembly = sp[11]
                if assembly == 'GRCh37':
                    if rs in rs2ss:
                        if ss != rs2ss[rs][0]:
                            problems[rs] = True
                            #print 'problem', ss, rs, rs2ss[rs]
                    else:
                        rs2ss[rs] = (ss, flip_bases)
    for rs in problems:
        del rs2ss[rs]
    return rs2ss
                        
def main():
    """Entry point"""

    coverage_cut = 8
    quality_cut = float(50)
    # 1 means flip
    rs2ss_flip = check_dbsnp_illumina()
    snps = load_exome_snps(quality_cut, coverage_cut)
    cmp_snpchip(snps['yusik'],
                snps['yusikPLX'],
                rs2ss_flip)    
    
if __name__ == '__main__':
    main()

