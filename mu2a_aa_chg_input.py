"""Make input for mu2a by ignoring yuaker
   and pooling aa_chg"""
import global_settings, mutations, os

def get_mutations_for_paired_samples(quality_cutoff, coverage_cutoff):
    """Load mutations from data/all_non_ref_hg19 for all paired samples"""

    sample2mutations = {}
    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg19/'

    for cancer, normal in global_settings.pairs:
        sample_name = cancer.split('0')[0]
        exome2mutations = {}
    
        cancer_qualities = mutations.get_consensus_qualities(os.path.join(use_data_dir, cancer, 
                                                                          sample_name + 'T.ann'))
        normal_qualities = mutations.get_consensus_qualities(os.path.join(use_data_dir, cancer,
                                                                          sample_name + 'N.ann'))

        for exome in global_settings.exome_types:
            data_file = os.path.join(use_data_dir, cancer, exome)
            inherited, somatic, murim = mutations.get_mutations(data_file, normal_qualities,
                                                                cancer_qualities, quality_cutoff,
                                                                False, coverage_cutoff)
            exome2mutations[exome] = (inherited, somatic)
        sample2mutations[sample_name] = exome2mutations

    return sample2mutations 

def write(afile, muts):
    """Write mutations to working/file"""

    with open(afile, 'w') as f:
        for chrpos in muts:
            chr, pos = chrpos.split(':')
            f.write(chr + '\t' + pos + '\n')

def screen_mutation_types(mutation_ls):
    """Allow mutations seen in the somatic mutation melanoma paper"""

    allowed = {}
    for chrpos in mutation_ls:
         mutation_type = mutation_ls[chrpos][0]
         if mutation_type not in ('AB:AA', 'AB:BB'): #AA:BB AA:AB?
              allowed[chrpos] = True
    return allowed

def main():
    """Call methods for somatic & inherited tables."""

    quality_cutoff = float(0) # really gt 50
    coverage_cutoff = int(8)

    sample2mutations = get_mutations_for_paired_samples(quality_cutoff, 
                                                        coverage_cutoff)
    inherited = {}
    somatic = {}
    for sample in sample2mutations:
        if 'yuaker' not in sample:
            i, s = sample2mutations[sample]['exome.aa_chg']
            sample_inherited = screen_mutation_types(i[sample])
            sample_somatic = screen_mutation_types(s[sample])
            for chrpos in sample_inherited:
                inherited[chrpos] = True
            for chrpos in sample_somatic:
                somatic[chrpos] = True
    write('working/inherited.mu2a', inherited)
    write('working/somatic.mu2a', somatic)
    print len(set(somatic.keys()) & set(inherited.keys()))

if __name__ == '__main__':
    main()
