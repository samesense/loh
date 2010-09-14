"""I'm trying to figure out why we have ~200 more mutations that Murim"""
import os, sys, global_settings
import cmp_murim_mutations_yusan, mutations

def get_my_mutations(quality_cutoff, coverage_cutoff):
    """Load mutations from working/"""

    # my_mutations = {}
    # with open('/home/perry/Projects/loh/working/murim.exome.aa_chg.vars') as f:
    #     for line in f:
    #         my_mutations[line.strip()] = True
    # return my_mutations

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
        for s in somatic['yusan']: all_somatic[s] = True
        for i in inherited['yusan']: all_inherited[i] = True
    return (all_somatic, all_inherited)

def get_murim_mutations(quality_cutoff):
    """Get somatic (not coding) mutations above some quality"""

    os.system("grep somatic /home/perry/Projects/loh/data/murim/CANCER_specific.csv | cut -f 3,4,10 > /home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual")
    murim_mutations = {}
    with open('/home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual') as f:
        for line in f:
            sp = line.strip().split('\t')
            if float(sp[-1]) > float(quality_cutoff): # murim uses 80
                murim_mutations[sp[0] + ':' + sp[1]] = True
    return murim_mutations

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

if __name__ == '__main__':
    main()
    
