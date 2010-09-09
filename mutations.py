"""Use the all_non_ref files to find somatic and inherited mutations"""
import global_settings, os
from collections import defaultdict

def init_zero(): return 0

def get_mutations(afile):
    """Count mutations AA:BB, AB:AA ... from raw data file"""

    inherited = defaultdict(init_zero)
    somatic = defaultdict(init_zero)
    with open(afile) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[9].split('-')
            idx = 18
            chr = sp[2]
            pos = sp[3]
            mutation_type = sp[4]
            normal, cancer = samples
            sample_name_normal = normal[0:-1]
            sample_name_cancer = cancer[0:-1]
            if sample_name_normal != sample_name_cancer:
                i = 1/0
            (avg_quality, cancer_quality, 
             normal_quality) = (float(x) for x in sp[14:17])
            if cancer_quality > float(100) and normal_quality > float(100):
                if mutation_type == 'AA:AA':
                    pass
                elif mutation_type in ('BB:BB', 'AB:AB'):
                    inherited[sample_name_normal] += 1
                else:
                    somatic[sample_name_normal] += 1
    return (inherited, somatic)

for exome_type in global_settings.exome_types:
    data_file = os.path.join('data/all_non_ref/',
                             exome_type)
    inherited, somatic = get_mutations(data_file)
    for sample in inherited:
        i = inherited[sample]
        s = somatic[sample]
        print('%s\t%s\t%d\t%d\t%.2f' %
              (exome_type, sample, i, s, 
               float(100)*float(s)/float(i+s)))
            
