"""I'm trying to figure out why we have ~200 more mutations that Murim"""
import os, sys, global_settings
import cmp_murim_mutations_yusan, mutations, bed_tools

def get_murim_covered(quality_cut):
    """Load chr locations from EX_CANC_1ln.snpfilter_anno_shortall_v1.62.txt for cancer and normal. If he doesn't have calls for both cancer and normal, I'm throwing it out of the comparison"""

    murim_data_dir = 'data/murim/'
    murim_suffix_normal = 'EX_BLD1_1ln.snpfilter_anno_shortall_v1.62.txt'
    murim_suffix_cancer = 'EX_CANC_1ln.snpfilter_anno_shortall_v1.62.txt'
    murim_normal = cmp_murim_mutations_yusan.load_murims_calls(os.path.join(murim_data_dir,
                                                  'yusanN',
                                                  murim_suffix_normal), 
                                     quality_cut)
    murim_cancer = cmp_murim_mutations_yusan.load_murims_calls(os.path.join(murim_data_dir,
                                                  'yusanT',
                                                  murim_suffix_cancer), 
                                     quality_cut)
    return set(murim_cancer.keys()) & set(murim_normal.keys())    

def get_my_mutations(quality_cutoff, coverage_cutoff):
    """Load mutations from working/"""

    # my_mutations = {}
    # with open('/home/perry/Projects/loh/working/murim.exome.aa_chg.vars') as f:
    #     for line in f:
    #         my_mutations[line.strip()] = True
    # return my_mutations

    bed_file = 'data/nimblegen/2.1M_Human_Exome_Annotation/2.1M_Human_Exome.bed'
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
        # only use the bed_tools NimbleGen
        # restriction for hg18 data
        for s in somatic['yusan']: 
            chr, pos = s.split(':')
            if bed_tools.find_location_in_bed(chr, int(pos), 
                                              bed_chr2posLs,
                                              bed_chr2st2end):
                all_somatic[s] = True
        for i in inherited['yusan']: 
            chr, pos = s.split(':')
            if bed_tools.find_location_in_bed(chr, int(pos), 
                                              bed_chr2posLs,
                                              bed_chr2st2end):
                all_inherited[i] = True
    return (set(all_somatic.keys()) & set(get_murim_covered(quality_cutoff)), set(all_inherited.keys()) & set(get_murim_covered(quality_cutoff)))

def remove_deletions(murim_mutations):
    """Trying to figure out why we are different from Murim.
       I see that we don't have calls for positions because
       the calls are bad; mostly deletions.
       I'm removing these from Murim's mutations"""

    pass

def get_murim_mutations(quality_cutoff):
    """Get somatic (not coding) mutations above some quality"""

    os.system("grep somatic /home/perry/Projects/loh/data/murim/CANCER_specific.csv | cut -f 3,4,10 > /home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual")
    
    murim_mutations = {}
    with open('/home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual') as f:
        for line in f:
            sp = line.strip().split('\t')
            if float(sp[-1]) > float(quality_cutoff): # murim uses 80
                murim_mutations[sp[0] + ':' + sp[1]] = True
    remove_deletions(murim_mutations)
    return set(murim_mutations.keys()) & set(get_murim_covered(quality_cutoff))

#print 'both', len(set(my_mutations) & set(murim_mutations))
#print 'us', len(my_mutations)
#print 'murim', len(murim_mutations)

def get_mutations_we_miss(quality_cutoff):
    """Return chrpos not in our results"""

    murim_mutations = get_murim_mutations(quality_cutoff)
    my_mutations = get_my_mutations(quality_cutoff)

    missing = {}
    for chrpos in murim_mutations:
        if chrpos not in my_mutations:
            missing[chrpos] = True
    return missing

def get_our_extra_mutations(quality_cutoff):
    """Return chrpos not Murim's results"""

    murim_mutations = get_murim_mutations(quality_cutoff)
    my_mutations = get_my_mutations(quality_cutoff)

    missing = {}
    for chrpos in my_mutations:
        if chrpos not in murim_mutations:
            missing[chrpos] = True
    return missing

def main():
    """Script entry"""

    #for chrpos in get_our_extra_mutations(float(sys.argv[1])):
    #    print chrpos

    # how many of 565 (no quality cutoff) do we get?
    quality = float(sys.argv[1])
    coverage = int(sys.argv[2])
    murim_mutations = set(get_murim_mutations(quality))
    our_somatic_mutations, our_inherited_mutations = [set(x) for x in get_my_mutations(quality, coverage)]
    print 'murim somatic\tus somatic\tus inherited\tmurim & us somatic\tmurim somatic & us inherited'
    print len(murim_mutations), len(our_somatic_mutations), len(our_inherited_mutations), len(murim_mutations & our_somatic_mutations),  len(murim_mutations & our_inherited_mutations)

    with open('working/yong_missing_from_murim', 'w') as f:
        for chrpos in murim_mutations:
            if chrpos not in our_somatic_mutations:
                f.write(chrpos + '\n')
    with open('working/yong_extra', 'w') as f:
        for chrpos in our_somatic_mutations:
            if chrpos not in murim_mutations:
                f.write(chrpos + '\n')

if __name__ == '__main__':
    main()
    
