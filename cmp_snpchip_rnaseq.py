"""Look at recall of SNP chip by RNA-seq"""
import cmp_snpchip_exomeseq, mu2a

def eval_rnaseq(patient, snpchip_calls):
    """Compare RNA seq call to snp chip call"""

    # RNA seq mutations that correspond to SNPs
    mu2a_file = 'data/rna_seq/rna_mutations_Mu2aOutput.txt'
    mu2a_output = mu2a.MU2A(mu2a_file)
    mu2a_snps = mu2a_output.get_snp_positions()
    print len(mu2a_snps.keys())
    gyang_file = 'data/rna_seq/rna_result.csv'
    with open(gyang_file) as f:
        for line in f:
            if patient in line.lower():
                sp = line.split('\t')
                chr = sp[7]
                pos = sp[8]
                chrpos = chr + ':' + pos
                if chrpos in mu2a_snps:
                    print chrpos
                #else:
                #    print 'missing', chrpos
def main():
    """Entry point"""

    snp_chip_file_cancer = 'working/yurif_cancer.call'
    snpchip_cancer = cmp_snpchip_exomeseq.get_snpchip_calls(snp_chip_file_cancer)
    eval_rnaseq('yusik', snpchip_cancer)

if __name__ == '__main__':
    main()
