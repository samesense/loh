"""Use the all_non_ref files to find somatic and inherited mutations"""
import global_settings, os
from collections import defaultdict

def init_zero(): return 0

def get_mutations(afile):
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
            if sample_name_normal != sample_name_cancer:
                i = 1/0
            (avg_quality, cancer_quality, 
             normal_quality) = (float(x) for x in sp[14:17])
            if cancer_quality > float(100) and normal_quality > float(100) and normal_coverage >= 8 and cancer_coverage >= 8:
                if mutation_type == 'AA:AA':
                    pass
                elif mutation_type in ('BB:BB', 'AB:AB'):
                    inherited[sample_name_normal][chrpos] = True
                else:
                    somatic[sample_name_normal][chrpos] = True
                    # murim can only see mutations 
                    # that differ from the reference
                    # only the first case matters
                    # at this quality cutoff
                    if mutation_type not in ('AB:AA',  'AA:AB', 'AA:BB'):
                        murim[sample_name_normal][chrpos] = True

    return (inherited, somatic, murim)

total_somatic_mutations = defaultdict(init_zero)
total_murim_mutations = defaultdict(init_zero)
for exome_type in global_settings.exome_types:
    data_file = os.path.join('data/all_non_ref/',
                             exome_type)
    inherited, somatic, murim = get_mutations(data_file)
    for sample in inherited:
        i = len(inherited[sample].keys())
        s = len(somatic[sample].keys())
        m = len(murim[sample].keys())
        total_somatic_mutations[sample] += s
        total_murim_mutations[sample] += m
        print('%s\t%s\t%d\t%d\t%d\t%.2f' %
              (exome_type, sample, i, s, m,
               float(100)*float(s)/float(i+s)))
for sample in total_somatic_mutations:
    print 'total somatic', sample, total_somatic_mutations[sample]
    print 'total murim', sample, total_murim_mutations[sample]
