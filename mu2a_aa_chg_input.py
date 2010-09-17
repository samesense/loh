"""Make input for mu2a by ignoring yuaker
   and pooling aa_chg"""
import call_class

# def get_mutations_for_paired_samples(quality_cutoff, coverage_cutoff):
#     """Load mutations from data/all_non_ref_hg19 for all paired samples"""

#     sample2mutations = {}
#     use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg19/'

#     for cancer, normal in global_settings.pairs:
#         sample_name = cancer.split('0')[0]
#         exome2mutations = {}
    
#         cancer_qualities = mutations.get_consensus_qualities(os.path.join(use_data_dir, cancer, 
#                                                                           sample_name + 'T.ann'))
#         normal_qualities = mutations.get_consensus_qualities(os.path.join(use_data_dir, cancer,
#                                                                           sample_name + 'N.ann'))

#         for exome in global_settings.exome_types:
#             data_file = os.path.join(use_data_dir, cancer, exome)
#             inherited, somatic, murim = mutations.get_mutations(data_file, normal_qualities,
#                                                                 cancer_qualities, quality_cutoff,
#                                                                 False, coverage_cutoff)
#             exome2mutations[exome] = (inherited, somatic)
#         sample2mutations[sample_name] = exome2mutations

#     return sample2mutations 

# def write(afile, muts, file_type):
#     """Write mutations to working/file"""

#     with open(afile, 'w') as f:
#         for chrpos in muts:
#             chr, pos = chrpos.split(':')
#             normal_call, cancer_call, ref_call = muts[chrpos]
#             if file_type == 'inherited':
#                 f.write(chr + '\t' + pos + '\t' 
#                         + ref_call + '\t' + cancer_call + '\n')
#             elif file_type == 'somatic':
#                 f.write(chr + '\t' + pos + '\t' 
#                         + normal_call + '\t' + cancer_call + '\n')
#             else:
#                 raise ValueError

# def screen_mutation_types(mutation_ls):
#     """Allow mutations seen in the somatic mutation melanoma paper"""

#     allowed = {}
#     for chrpos in mutation_ls:
#          mutation_type, normal_call, cancer_call, ref_allele = mutation_ls[chrpos]
#          if mutation_type not in ('AB:AA', 'AB:BB'): #AA:BB AA:AB?
#               allowed[chrpos] = (normal_call, cancer_call, ref_allele)
#     return allowed

def main():
    """Call methods for somatic & inherited tables."""

    quality_cutoff = float(0) # really gt 50
    coverage_cutoff = int(8)

    sample2exome2mutations = call_class.get_mutations_for_paired_samples(quality_cutoff,
                                                                         coverage_cutoff)

    # accumulate across samples
    inherited = {}
    somatic = {}
    for sample in sample2exome2mutations:
        if 'yuaker' not in sample:
            sample_inherited, sample_somatic = sample2exome2mutations[sample]
            for chrpos in sample_inherited['exome.aa_chg']:
                inherited[chrpos] = sample_inherited['exome.aa_chg'][chrpos]
            for chrpos in sample_somatic['exome.aa_chg']:
                somatic[chrpos] = sample_somatic['exome.aa_chg'][chrpos]
    call_class.write_MU2A_input(inherited, 'working/mu2a/inherited.mu2a_input.new', 'inherited')
    call_class.write_MU2A_input(somatic, 'working/mu2a/somatic.mu2a_input.new', 'somatic')

    # inherited = {}
    # somatic = {}
    # for sample in sample2mutations:
    #     if 'yuaker' not in sample:
    #         i, s = sample2mutations[sample]['exome.aa_chg']
    #         sample_inherited = screen_mutation_types(i[sample])
    #         sample_somatic = screen_mutation_types(s[sample])
    #         for chrpos in sample_inherited:
    #             inherited[chrpos] = sample_inherited[chrpos]
    #         for chrpos in sample_somatic:
    #             somatic[chrpos] = sample_somatic[chrpos]
    # write('working/mu2a/inherited.mu2a_input', inherited, 'inherited')
    # write('working/mu2a/somatic.mu2a_input', somatic, 'somatic')
    # print len(set(somatic.keys()) & set(inherited.keys()))

if __name__ == '__main__':
    main()
