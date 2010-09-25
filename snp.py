"""Tools for getting RS SNPs from NCBI,
   and figuring out nucleutides from
   Illumina rs# AB calls"""
import os, random, sys, global_settings

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

    with open('working/snps', 'w') as outf:
        with open('data/dbsnp/seq_snp.md') as f:
            for line in f:
                (tax_id, chromosome, chr_start,
                 chr_stop, chr_orient, contig,  
                 ctg_start, ctg_stop, ctg_orient,
                 feature_name, feature_id, 
                 feature_type, group_label, 
                 weight) = line.strip().split('\t')
                if feature_type == 'SNP' and group_label == 'Primary_Assembly' and feature_name in snp_tb:
                    outf.write('%s\t%s\t%s\t%s\t%s\t%s\n' %
                               (feature_name, chromosome, chr_start,
                                chr_stop, chr_orient, snp_tb[feature_name]))

def get_TB_for_SNPs(use_ss):
    """For all ss SNPs, get TOP|BOTTOM"""

    snps = {}
    with open('data/dbsnp/SubSNP.bcp') as f:
        for line in f:
            sp = line.strip().split('\t')
            ss_snp = 'ss' + sp[0]
            if ss_snp in use_ss:
                tb = sp[-4]
                snps[ss_snp] = tb
                #print tb, sp[-4], sp[-6]
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
                alleles = sp[2]
                rs = sp[4]
                ss2rs_orien = sp[5]
                rs2genome_orien = sp[10]
                assembly = sp[11]
                if assembly == 'GRCh37':
                    if rs in rs2ss:
                        if ss != rs2ss[rs][0]:
                            problems[rs] = True
                            #print 'problem', ss, rs, rs2ss[rs]
                    else:
                        flip_bases = False
                        if ss2rs_orien == '1' or rs2genome_orien == '1':
                            if not (ss2rs_orien == '1' and 
                                    rs2genome_orien == '1'):
                                flip_bases = True
                        rs2ss[rs] = [ss, alleles, flip_bases]
    ss = {}
    for rs in problems:
        del rs2ss[rs]
    for rs in rs2ss:
        ss[rs2ss[rs][0]] = True
    return (rs2ss, ss)

def check_TB_problems(rs2ss, ss2tb):
    """See if there is disagreement between alleles and top|bottom"""

    #problems = {}
    #total = {}
    for rs in rs2ss:
        ss, alleles, flip_bases = rs2ss[rs]
        if ss in ss2tb:
            #total[ss] = True
            # append TB
            rs2ss[rs].append(ss2tb[ss])
            if alleles == 'A/G' or alleles == 'G/A' or alleles == 'A/C' or alleles == 'C/A':
                if ss2tb[ss] != 'T':
                    print ss, alleles, ss2tb[ss]
            elif alleles == 'T/C' or alleles == 'C/T' or alleles == 'T/G' or alleles == 'G/T':
                if ss2tb[ss] != 'B':
                    print 'here', ss, alleles, ss2tb[ss]
#    print 'problems', len(problems), len(total)

# def get_TB_for_SNPs():
#     """For all rs SNPs get TOP|BOTTOM"""
    
#     snps = {}
#     with open('data/dbsnp/SubSNP_top_or_bot.bcp') as f:
#         for line in f:
#             (snp, tb, stuff, datetime) = line.split('\t')
#             snps['rs' + snp] = tb
#     return snps
                
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

def get_genotype_for_SNPs():
    """Open working/snps and append genotype"""

    snps = {}
    with open('working/snps') as oldf:
        for line in oldf:
            (snp, chr, st, stp, orient, tb) = line.strip().split('\t')
            snps[snp] = (chr, st, stp, orient, tb)
    chrs =  [str(x) for x in range(1,22)]
    chrs.extend(['X', 'Y'])
    with open('working/snps_full', 'w') as outf:
        for chr in chrs:
            with open('data/dbsnp/chr%s' % chr) as f:
                for line in f:
                    sp = line.split('|')
                    snp = sp[2].split()[0]
                    genotype = sp[8].split('=')[1].strip("\"")
                    if snp in snps:
                        (chr, st, stp, orient, tb) = snps[snp]
                        outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                   (snp, chr, st, stp, orient, tb, genotype))

def flip_if_needed(alleleA, alleleB, flip):
    """Take the complement if needed"""

    if flip:
        return (global_settings.comp[alleleA], 
                global_settings.comp[alleleB])
    else:
        return (alleleA, alleleB)

def get_allele_for_genotype(snp_info, GType, SNP):
    """Given Illumina's call, find the nucleutodes"""

    if GType == 'NC':
        sys.stderr.write('err2 ' + GType + '\n')
        return ''
    else:
        alleleA = ''
        alleleB = ''
        (ss, alleles, flip_bases, tb) = snp_info
        if tb == 'T':
            if alleles == 'A/G' or alleles == 'G/A':
                alleleA = 'A'
                alleleB = 'G'
            elif alleles == 'A/C' or alleles == 'C/A':
                alleleA = 'A'
                alleleB = 'C'
            elif alleles == 'A/T' or alleles == 'T/A':
                alleleA = 'A'
                alleleB = 'T'
            elif alleles == 'C/G' or alleles == 'G/C':
                alleleA = 'C'
                alleleB = 'G'
            else:
                sys.stderr.write('err1 ' + SNP + ' ' + str(snp_info) + '\n')
                #raise ValueError
        elif tb == 'B':
            if alleles == 'T/C' or alleles == 'C/T':
                alleleA = 'T'
                alleleB = 'C'
            elif alleles == 'T/G' or alleles == 'G/T':
                alleleA = 'T'
                alleleB = 'G'
            elif alleles == 'A/T' or alleles == 'T/A':
                alleleA = 'T'
                alleleB = 'A'
            elif alleles == 'C/G' or alleles == 'G/C':
                alleleA = 'G'
                alleleB = 'C'
            else:
                sys.stderr.write('err1 ' + SNP + ' ' + str(snp_info) + '\n')
                #raise ValueError    
        else:
            # why is t/b missing???
            sys.stderr.write('err1 ' + SNP + ' ' + str(snp_info) + '\n')
            #raise ValueError

        if alleleB and alleleB:
            if GType == 'AA':
                return flip_if_needed(alleleA, alleleA, flip_bases)
            elif GType == 'AB':
                return flip_if_needed(alleleA, alleleB, flip_bases)
            elif GType == 'BB':
                return flip_if_needed(alleleB, alleleB, flip_bases)
            else:
                sys.stderr.write('err2 ' + GType + '\n')
                return ''
                #raise ValueError
    
def convert_snpchip(snpchip_file, rs2ss_info):
    """Take Illumina's AB calls and convert to nucleotides"""

    # snp_info = {}
    # with open('working/snps_full') as f:
    #     for line in f:
    #         (snp, chr, st, stp, 
    #          orient, tb, alleles) = line.strip().split('\t')
    #         snp_info[snp] = (chr, st, stp, 
    #                          orient, tb, alleles)
    with open(snpchip_file) as f:
        f.readline()
        for line in f:
            (Index, Name, Address, Chr, Position, GenTrain,
             A, C, G, T, GType, Score, Theta, R) = line.strip().split('\t')
            if Name in rs2ss_info:
                call = get_allele_for_genotype(rs2ss_info[Name],
                                               GType, Name)
                if call:
                    print('%s\t%s\t%s' % 
                          (Name, GenTrain, '/'.join(call)))

def main():
    """Entry point"""

    rs2ss, ss = check_dbsnp_illumina()
    ss2tb = get_TB_for_SNPs(ss)
    check_TB_problems(rs2ss, ss2tb)
    convert_snpchip('data/snp_chip/yuiri_tumor',
                    rs2ss)

    # get_TB_strand_for_SNPs()
    # get_genotype_for_SNPs()
    # convert_snpchip('data/snp_chip/yuiri_normal')

if __name__ == '__main__':
    main()

