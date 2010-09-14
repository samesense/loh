"""Find out why we don't recover murim's mutations"""
import sys, nose.tools
sys.path.append('../')
import mutations, cmp_murim_mutations_yusan_2

def test_murim_recall_annFiles():
    """See if we have all Murim's mutations as somatic"""

    we_are_missing = cmp_murim_mutations_yusan_2.get_mutations_we_miss()
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
 

    
