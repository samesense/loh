"""Tools for getting RS SNPs from NCBI,
   and figuring out nucleutides from
   Illumina rs# AB calls"""
import os, random

random.seed()

def get_rs_snps(snp_ls):
    """Get RS# SNP data from NCBI"""

    tmp_snp = 'tmpSNP' + str(random.randint(0,1000))
    query = "'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=" + ','.join(snp_ls) + "&rettype=flt&email=samesense@gmail.com'"
    print query
    os.system('wget -O ' + tmp_snp +
              ' ' + query)

def get_orientation_SNPs(snp_tb):
    """For all rs SNPs, get orientation (+/-) and location"""

    with open('data/dbsnp/seq_snp.md') as f:
        for line in f:
            (tax_id, chromosome, chr_start,
             chr_stop, chr_orient, contig,  
             ctg_start, ctg_stop, ctg_orient,
             feature_name, feature_id, 
             feature_type, group_label, 
             weight) = line.strip().split('\t')
            if feature_type == 'SNP' and group_label == 'HuRef' and feature_name in snp_tb:
                print('%s\t%s\t%s\t%s\t%s\t%s' %
                      (feature_name, chromosome, chr_start,
                       chr_stop, chr_orient, snp_tb[feature_name]))

def get_TB_for_SNPs():
    """For all rs SNPs get TOP|BOTTOM"""
    
    snps = {}
    with open('data/dbsnp/SubSNP_top_or_bot.bcp') as f:
        for line in f:
            (snp, tb, stuff, datetime) = line.split('\t')
            snps['rs' + snp] = tb
    return snps
                
def get_TB_strand_for_SNPs():
    """For all rs SNPs, get TOP|BOTTOM and +/'"""

    snp_tb = get_TB_for_SNPs()
    snp_orientation = get_orientation_SNPs(snp_tb)

    # for snp in snp_orientation:
    #     if snp in snp_tb:
    #         snp_orientation[snp].append(snp_tb[snp])
    # print len(snp_orientation)
    # for snp in snp_orientation:
    #     if snp not in snp_tb:
    #         del snp_orientation[snp]

    # print len(snp_orientation)

def main():
    get_TB_strand_for_SNPs()

if __name__ == '__main__':
    main()

