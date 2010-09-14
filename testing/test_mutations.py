"""Find out why we don't recover murim's mutations"""
import sys, nose.tools
sys.path.append('../')
import mutations, cmp_murim_mutations_yusan_2

def test_murim_recall_annFiles():
    """See if we have all Murim's mutations as somatic"""

    quality_cutoff = float(100)
    we_are_missing = cmp_murim_mutations_yusan_2.get_mutations_we_miss(quality_cutoff)
    use_data_dir = '../data/all_non_ref_hg18/'
    cancer_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanT.ann')
    normal_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanN.ann')
    cancer_intr = len(we_are_missing) - len(set(cancer_qualities) & set(we_are_missing))
    normal_intr = len(we_are_missing) - len(set(normal_qualities) & set(we_are_missing))
    for chrpos in we_are_missing:
        nose.tools.assert_true(chrpos in normal_qualities,
                               chrpos + ' missing in normal, ' + str(cancer_intr) + ' ' + str(normal_intr))
        nose.tools.assert_true(chrpos in cancer_qualities,
                               chrpos + ' missing in cancer, ' + str(cancer_intr) + ' ' + str(normal_intr))

def test_AA_AA_calls():
    """See if the AA:AB calls are screwed up. See if we see a change when there is none."""

    quality_cutoff = float(100)
    extra_mutations = cmp_murim_mutations_yusan_2.get_our_extra_mutations(quality_cutoff)
    use_data_dir = '../data/all_non_ref_hg18/'
    cancer_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanT.ann')
    normal_qualities = mutations.get_consensus_qualities(use_data_dir + 'yusanN.ann')
    data_file = '../data/all_non_ref_hg18/exome.aa_chg'
    inherited, somatic, murim = mutations.get_mutations(data_file, normal_qualities,
                                                        cancer_qualities, quality_cutoff,
                                                        True)
    for chrpos in extra_mutations:
        mutation_type, normal_call, cancer_call = somatic['yusan'][chrpos]
        m1, m2 = mutation_type.split(':')
        nose.tools.assert_true(m1 != m2,
                               mutation_type + ' should not be here')

# def test_murim_recall_somatic():
#     """See how many of Murim's 565 somatic mutations we recover"""

#     quality_cutoff = float(-1)
#     murim_mutations = get_mutations_
 

    
