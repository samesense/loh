"""
Table for yusan

inherited

-Number of HQ SNVs
-Number of novel HQ SNVs (ie. not in dbSNP)
-Number of novel HQ homoz SNVs
-Number of novel HQ heteroz SNVs.

somatic

-number of high quality somatic mutations
-Number of novel HQ homoz somatic mutations
-Number of novel HQ heteroz somatic mutations.
-How many of the somatic mutations are AA changing
-How many are affecting domains
-maybe some other stats...
"""
import sys, os
import mutations, global_settings, bed_tools, call_class

def get_mutations_yusan(quality_cutoff, coverage_cutoff):
    """Load mutations from data/all_non_ref_hg18"""


    exome2mutations = {}
    bed_file = 'data/nimblegen/2.1M_Human_Exome_Annotation/2.1M_Human_Exome.bed'
    bed_chr2st2end, bed_chr2posLs = bed_tools.load_bed(bed_file, 
                                                       'NimbleGen Tiled Regions')
    # NimbleGen Tiled Regions
    # Target Regions

    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg18/'
    all_somatic = {}
    all_inherited = {}
    cancer_qualities = call_class.get_consensus_qualities(use_data_dir + 'yusanT.ann')
    normal_qualities = call_class.get_consensus_qualities(use_data_dir + 'yusanN.ann')
    for exome in global_settings.exome_types:
        data_file = use_data_dir + exome
        inherited, somatic, murim = mutations.get_mutations(data_file, normal_qualities,
                                                            cancer_qualities, quality_cutoff,
                                                            False, coverage_cutoff)
        exome2mutations[exome] = (inherited, somatic)
    return exome2mutations       

def screen_mutation_types(mutation_ls):
    """Allow mutations seen in the somatic mutation melanoma paper"""

    allowed = {}
    for chrpos in mutation_ls:
         mutation_type = mutation_ls[chrpos][0]
         if mutation_type not in ('AB:AA', 'AB:BB'): #AA:BB AA:AB?
              allowed[chrpos] = True
    return allowed

def intersection_table(exome2mutations):
    """Print intersection between aa_chg normal->cancer mutations"""
    
    exome = 'exome.aa_chg'
    print 'Sample\tExome\tInherited count\tSomatic count'
    for sample1 in exome2mutations:
        for sample2 in exome2mutations:
            if sample1 != sample2 and sample1 != 'yuaker' and sample2 != 'yuaker':
                inherited_pre_s1, somatic_pre_s1 = exome2mutations[sample1][exome]
                somatic1 = screen_mutation_types(somatic_pre_s1[sample1])

                inherited_pre_s2, somatic_pre_s2 = exome2mutations[sample2][exome]
                somatic2 = screen_mutation_types(somatic_pre_s2[sample2])

                both = len(set(somatic2) & set(somatic1))
                print('%s\t%s\t%d\t%d\t%d' %
                      (sample1, sample2, 
                       len(somatic1.keys()),
                       len(somatic2.keys()),
                       both))

def somatic_table_melanoma_paper_paired_samples(sample2mutations):
    """
    Discard cases where the tumor allele is present in normal,
    or when the tumor allele is the same as the reference
    """
    
    print 'Sample\tExome\tInherited count\tSomatic count'
    for sample in sample2mutations:
        inherited, somatic = sample2mutations[sample]
        for exome in set(inherited.keys()) | set(somatic.keys()):
            print('%s\t%s\t%d\t%d' %
                  (sample, exome, len(inherited[exome].keys()),
                   len(somatic[exome].keys()))) 

def somatic_table_melanoma_paper(exome2mutations):
    """
    Discard cases where the tumor allele is present in normal,
    or when the tumor allele is the same as the reference
    """
    
    print 'Exome\tInherited count\tSomatic count'
    for exome in exome2mutations:
        inherited_pre, somatic_pre = exome2mutations[exome]
        somatic = screen_mutation_types(somatic_pre['yusan'])
        inherited = screen_mutation_types(inherited_pre['yusan'])
        print('%s\t%s\t%s' %
              (exome, len(inherited.keys()),
               len(somatic.keys())))    

def somatic_table(exome2mutations):
    """
    -number of high quality somatic mutations
    -Number of novel HQ homoz somatic mutations
    -Number of novel HQ heteroz somatic mutations.
    -How many of the somatic mutations are AA changing
    -How many are affecting domains
    """
    
    print 'Exome\tInherited count\tSomatic count'
    for exome in exome2mutations:
        inherited, somatic = exome2mutations[exome]
        print('%s\t%d\t%d' %
              (exome, len(inherited['yusan'].keys()),
               len(somatic['yusan'].keys())))

def main():
    """Call methods for somatic & inherited tables."""

    quality_cutoff = float(0)
    coverage_cutoff = int(8)

    sample2exome2mutations = call_class.get_mutations_for_paired_samples(quality_cutoff,
                                                                         coverage_cutoff)
    somatic_table_melanoma_paper_paired_samples(sample2exome2mutations)

#    intersection_table(exome2mutations)

#    exome2mutations = get_mutations_yusan(quality_cutoff, coverage_cutoff)
#    somatic_table_melanoma_paper(exome2mutations)

    

    #mu2a_somatic = mu2a.MU2A('working/mu2a/somatic.mu2a.output')
    #print len(mu2a_somatic.get_kinse_mutations)
    
#somatic_table(exome2mutations)
    #print '++++++++++++++++++++'
    #

if __name__ == '__main__':
    main()
