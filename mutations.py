"""Use the all_non_ref files to find somatic and inherited mutations. Use consensus quality from *ann files to find quality cutoffs."""
import global_settings, os, sys
from collections import defaultdict

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
 
def get_mutations(afile, normal_qualities, cancer_qualities, quality_cutoff):
    """Count mutations AA:BB, AB:AA ... from raw data file"""

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
            if normal_qualities[chrpos] > quality_cutoff and cancer_qualities[chrpos] > quality_cutoff and normal_coverage >= 8 and cancer_coverage >= 8:
                if mutation_type == 'AA:AA':
                    pass
                elif mutation_type in ('BB:BB', 'AB:AB'):
                    inherited[sample_name_normal][chrpos] = (mutation_type,
                                                             normal_call,
                                                             cancer_call)
                else:
                    somatic[sample_name_normal][chrpos] = (mutation_type,
                                                           normal_call,
                                                           cancer_call)
                    # murim can only see mutations 
                    # that differ from the reference
                    # only the first case matters
                    # at this quality cutoff
                    if mutation_type not in ('AB:AA', 'AA:AB', 'AA:BB', 'BB:AA'):
                        murim[sample_name_normal][chrpos] = (mutation_type,
                                                             normal_call,
                                                             cancer_call)

    return (inherited, somatic, murim)

quality_cutoff = float(100)

# grab consensus qualities from *ann files
cancer_qualities = get_consensus_qualities('data/all_non_ref/yusanT.ann')
normal_qualities = get_consensus_qualities('data/all_non_ref/yusanN.ann')

total_somatic_mutations = defaultdict(dict)
total_murim_mutations = defaultdict(dict)
total_somatic_exome_mutations = defaultdict(dict)
total_inherited_exome_mutations = defaultdict(dict)
for exome_type in global_settings.exome_types:
    data_file = os.path.join('data/all_non_ref/',
                             exome_type)
    inherited, somatic, murim = get_mutations(data_file, normal_qualities,
                                              cancer_qualities, quality_cutoff)
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

        print('%s\t%s\t%d\t%d\t%d\t%.2f' %
              (exome_type, sample, i, s, m,
               float(100)*float(s)/float(i+s)))
        if exome_type == 'exome.aa_chg':
            dump_mu2a_input(murim[sample])
sample = 'yusan'
dump_mutants(total_murim_mutations[sample], 'working/somatic.new_mutants')
dump_mutants(total_somatic_mutations[sample], 'working/total.new_mutants')
# no unknown or introns for these, which is why I call them coding
dump_mutants(total_somatic_exome_mutations[sample], 'working/somatic_coding.mutants')
dump_mutants(total_inherited_exome_mutations[sample], 'working/inherited_coding.mutants')

for sample in total_somatic_mutations:
    print 'total somatic', sample, len(total_somatic_mutations[sample])
    print 'total murim', sample, len(total_murim_mutations[sample])
    coding_somatic = len(total_somatic_exome_mutations[sample])
    coding_inherited = len(total_inherited_exome_mutations[sample])
    coding_total = float(coding_somatic + coding_inherited)
    print 'total somatic exome', sample, coding_somatic, float(100)*float(coding_somatic)/float(coding_total)
    print 'total inherited exome', sample, coding_inherited, float(100)*float(coding_inherited)/float(coding_total)

