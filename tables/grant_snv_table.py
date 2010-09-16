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
import sys
sys.path.append('../')
import mutations, global_settings, bed_tools

def get_mutations(quality_cutoff, coverage_cutoff):
    """Load mutations from data/all_non_ref_hg19"""


    exome2mutations = {}
    bed_file = '../data/nimblegen/2.1M_Human_Exome_Annotation/2.1M_Human_Exome.bed'
    bed_chr2st2end, bed_chr2posLs = bed_tools.load_bed(bed_file, 
                                                       'NimbleGen Tiled Regions')
    # NimbleGen Tiled Regions
    # Target Regions

    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg18/'
    all_somatic = {}
    all_inherited = {}
    cancer_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanT.ann')
    normal_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanN.ann')
    for exome in global_settings.exome_types:
        data_file = use_data_dir + exome
        inherited, somatic, murim = mutations.get_mutations(data_file, normal_qualities,
                                                            cancer_qualities, quality_cutoff,
                                                            False, coverage_cutoff)
        exome2mutations[exome] = (inherited, somatic)
    return exome2mutations        

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
        print('%s\t%s\t%s' %
              (exome, len(inherited['yusan'].keys()),
               len(somatic['yusan'].keys())))

def main():
    """Call methods for somatic & inherited tables."""

    quality_cutoff = float(100)
    coverage_cutoff = int(8)
    exome2mutations = get_mutations(quality_cutoff, coverage_cutoff)

    somatic_table(exome2mutations)

if __name__ == '__main__':
    main()
