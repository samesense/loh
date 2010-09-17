"""Make input for mu2a by ignoring yuaker
   and pooling aa_chg"""
import call_class

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
    call_class.write_MU2A_input(inherited, 'working/mu2a/inherited.mu2a_input', 'inherited')
    call_class.write_MU2A_input(somatic, 'working/mu2a/somatic.mu2a_input', 'somatic')

if __name__ == '__main__':
    main()
