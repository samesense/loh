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
    
