"""Use the all_non_ref files to find somatic and inherited mutations. Use consensus quality from *ann files to find quality cutoffs."""
import global_settings, os, sys
from collections import defaultdict
import cmp_murim_mutations_yusan

def get_consensus_qualities(afile):
    """Grab quality score for all calls. This is parsing the pileup output"""

    qualities = {}
    with open(afile) as f:
        for line in f:
            sp = line.strip().split('\t')
            chr = sp[2]
            pos = sp[3]
            con_quality = float(sp[13])
            qualities[chr+':'+pos] = con_quality
    return qualities

def dump_mu2a_input(chrpos2mutant):
    """Write input to mu2a"""

    with open('/home/perry/Dropbox/mu2a', 'w') as f:
        for chrpos in chrpos2mutant:
            chr, pos = chrpos.split(':')
            mutant, normal_call, cancer_call = chrpos2mutant[chrpos]
            f.write('%s\t%s\t%s\t%s\n' %
                    (chr, pos, normal_call, cancer_call))

def init_zero(): return 0

def dump_mutants(chrpos2mutant, dumpfile):
    """Write mutations to file"""

    with open(dumpfile, 'w') as f:
        for chrpos in chrpos2mutant:
            chr, pos = chrpos.split(':')
            mutant, normal_call, cancer_call = chrpos2mutant[chrpos]
            f.write('%s\t%s\t%s\t%s\t%s\n' %
                    (chr, pos, mutant, normal_call, cancer_call))
 
def get_mutations(afile, normal_qualities, cancer_qualities, quality_cutoff, cmp_murim, coverage_cutoff):
    """Count mutations AA:BB, AB:AA ... from raw data file"""

    if cmp_murim:
        limiting_locations_normal = cmp_murim_mutations_yusan.load_murims_calls('/home/perry/Projects/loh/data/murim/yusanN/EX_BLD1_1ln.snpfilter_anno_shortall_v1.62.txt', quality_cutoff)
        limiting_locations_cancer = cmp_murim_mutations_yusan.load_murims_calls('/home/perry/Projects/loh/data/murim/yusanT/EX_CANC_1ln.snpfilter_anno_shortall_v1.62.txt', quality_cutoff)
        limiting_locations = set(limiting_locations_normal.keys()) & set(limiting_locations_cancer.keys())

    inherited = defaultdict(dict)
    somatic = defaultdict(dict)
    murim = defaultdict(dict)

    with open(afile) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[9].split('-')
            idx = 18
            chr = sp[2]
            pos = sp[3]
            chrpos = chr + ':' + pos
            mutation_type = sp[4]
            normal, cancer = samples
            sample_name_normal = normal[0:-1]
            sample_name_cancer = cancer[0:-1]
            normal_coverage = int(sp[20])
            cancer_coverage = int(sp[25])
            normal_call = sp[18].upper()
            cancer_call = sp[23].upper()

            # throw exception if normal/cancer name does not match
            if sample_name_normal != sample_name_cancer:
                i = 1/0

            (avg_snp_quality, min_snp_quality, 
             max_snp_quality) = (float(x) for x in sp[14:17])
            if normal_qualities[chrpos] > quality_cutoff and cancer_qualities[chrpos] > quality_cutoff and normal_coverage >= coverage_cutoff and cancer_coverage >= coverage_cutoff and max_snp_quality > quality_cutoff:
                if cmp_murim and chrpos not in limiting_locations:
                    pass
                else:
                    if mutation_type == 'AA:AA':
                        pass
                    elif mutation_type in ('BB:BB', 'AB:AB'):
                        inherited[sample_name_normal][chrpos] = (mutation_type,
                                                                 normal_call,
                                                                 cancer_call)
                    elif normal_call != cancer_call:
                        somatic[sample_name_normal][chrpos] = (mutation_type,
                                                               normal_call,
                                                               cancer_call)
                        # murim can only see mutations 
                        # that differ from the reference
                        # only the first case matters
                        # at this quality cutoff
                        if mutation_type in ('AB:BB',):#not in ('AB:AA', 'AA:AB', 'AA:BB', 'BB:AA'): # these should give the same counts
                            murim[sample_name_normal][chrpos] = (mutation_type,
                                                                 normal_call,
                                                                 cancer_call)

    return (inherited, somatic, murim)

def main():
    """Call by default"""

    quality_cutoff = float(-1)
    coverage_cutoff = int(-1)
    use_data_dir = 'data/all_non_ref_hg18/' # | data/all_non_ref_hg19/
    cmp_murim = False
    # grab consensus qualities from *ann files
    cancer_qualities = get_consensus_qualities(use_data_dir + 'yusanT.ann')
    normal_qualities = get_consensus_qualities(use_data_dir + 'yusanN.ann')

    total_somatic_mutations = defaultdict(dict)
    total_inherited_mutations = defaultdict(dict)
    total_murim_mutations = defaultdict(dict)
    total_somatic_exome_mutations = defaultdict(dict)
    total_inherited_exome_mutations = defaultdict(dict)
    for exome_type in global_settings.exome_types:
        data_file = os.path.join(use_data_dir,
                                 exome_type)
        inherited, somatic, murim = get_mutations(data_file, normal_qualities,
                                                  cancer_qualities, quality_cutoff,
                                                  cmp_murim, coverage_cutoff)
        for sample in inherited:
            i = len(inherited[sample].keys())
            s = len(somatic[sample].keys())
            m = len(murim[sample].keys())

            if exome_type not in ('exome.intron', 'exome.unknown'):
                for chrpos in somatic[sample]:
                    total_somatic_exome_mutations[sample][chrpos] = somatic[sample][chrpos]
                for chrpos in inherited[sample]:
                    total_inherited_exome_mutations[sample][chrpos] = inherited[sample][chrpos]

            for chrpos in somatic[sample]:
                total_somatic_mutations[sample][chrpos] = somatic[sample][chrpos]

            for chrpos in murim[sample]:
                total_murim_mutations[sample][chrpos] = murim[sample][chrpos]

            for chrpos in inherited[sample]:
                total_inherited_mutations[sample][chrpos] = inherited[sample][chrpos]

            print('%s\t%s\t%d\t%d\t%d\t%.2f' %
                  (exome_type, sample, i, s, m,
                   float(100)*float(s)/float(i+s)))
            if exome_type == 'exome.aa_chg':
                dump_mu2a_input(murim[sample])
            with open('working/murim.' + exome_type + '.vars', 'w') as outf:
                for chrpos in murim[sample]:
                    outf.write(chrpos + '\n')
    sample = 'yusan'
    dump_mutants(total_murim_mutations[sample], 'working/somatic.new_mutants')
    dump_mutants(total_somatic_mutations[sample], 'working/total.new_mutants')
    # no unknown or introns for these, which is why I call them coding
    dump_mutants(total_somatic_exome_mutations[sample], 'working/somatic_coding.var')
    dump_mutants(total_inherited_exome_mutations[sample], 'working/inherited_coding.mutants')
    dump_mutants(total_somatic_mutations[sample], 'working/somatic_total.var')

    for sample in total_somatic_mutations:
        print 'total somatic', sample, len(total_somatic_mutations[sample])
        print 'total murim', sample, len(total_murim_mutations[sample])
        coding_somatic = len(total_somatic_exome_mutations[sample])
        coding_inherited = len(total_inherited_exome_mutations[sample])
        coding_total = float(coding_somatic + coding_inherited)
        total_somatic = len(total_somatic_mutations[sample])
        total_inherited = len(total_inherited_mutations[sample])
        all_total = total_somatic + total_inherited
        print 'total somatic exome', sample, coding_somatic, float(100)*float(coding_somatic)/float(coding_total)
        print 'total inherited exome', sample, coding_inherited, float(100)*float(coding_inherited)/float(coding_total)

        print 'total somatic', sample, total_somatic, float(100)*float(total_somatic)/float(all_total)
        print 'total inherited', sample, total_inherited, float(100)*float(total_inherited)/float(all_total)

if __name__ == '__main__':
    main()
