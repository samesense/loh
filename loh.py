"""Find somatic changes that are LOH in cancer."""
import global_settings, os, sys
from collections import defaultdict

def init_zero():
    return 0

def is_loh(cancer, normal):
    """Return true if normal->cancer is LOH"""

    if normal in global_settings.het_bases and cancer in global_settings.homo_bases:
        return True
    
    return False

def write_mutations(afile, of, q_limit):
    """Find instances where the call from normal to cancer changes in the initial exome files"""

    somatics = defaultdict(init_zero)
    totals = defaultdict(init_zero)
    with open(file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            calls = {}
            qualities = {}
            for s in samples:
                chr = sp[2].split('chr')[1]
                pos = sp[3]
                quality = float(sp[idx+1])
                coverage = int(sp[idx+2])
                call = sp[idx]
                calls[s] = call
                qualities[s] = quality
                idx += 5
            # check for changes in normal/cancer pairs
            for cancer, normal in global_settings.pairs:
                if cancer in calls and normal in calls:
                    if qualities[cancer] > float(q_limit) and qualities[normal] > q_limit:
                        totals[cancer+':'+normal] += 1
                        if calls[cancer] != calls[normal]:
                            if is_loh(calls[cancer], calls[normal]):
                                somatics[cancer+':'+normal] += 1
                                of.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                                         % (afile.split('/')[-1], chr, pos,
                                            normal, cancer,                                        
                                            calls[normal], calls[cancer]))
                    
    for sample in somatics:
        print afile.split('/')[-1] + '\t' + sample + '\t' + str(float(100)*float(somatics[sample])/float(totals[sample]))

def write_mutations_non_ref(afile, of, q_limit):
    """Find instances where the call from normal to cancer changes in the all_non_ref_hg19 files"""

    somatics = defaultdict(init_zero)
    totals = defaultdict(init_zero)
    with open(file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[9].split('-')
            idx = 18
            calls = {}
            qualities = {}
            for s in samples:
                chr = sp[2].split('chr')[1]
                pos = sp[3]
                quality = float(sp[idx+1])
                coverage = int(sp[idx+2])
                call = sp[idx]
                calls[s] = call
                qualities[s] = quality
                idx += 5
            # check for changes in normal/cancer pairs
            for cancer, normal in global_settings.pairs:
                if cancer in calls and normal in calls:
                    if qualities[cancer] > float(q_limit) and qualities[normal] > q_limit:
                        totals[cancer+':'+normal] += 1
                        if calls[cancer] != calls[normal]:
                            if is_loh(calls[cancer], calls[normal]):
                                somatics[cancer+':'+normal] += 1
                                of.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                                         % (afile.split('/')[-1], chr, pos,
                                            normal, cancer,                                        
                                            calls[normal], calls[cancer]))
                    
    for sample in somatics:
        print afile.split('/')[-1] + '\t' + sample + '\t' + str(float(100)*float(somatics[sample])/float(totals[sample]))

def main():
    """Entry point. Use both old data (exome) and new data(all_non_ref_hg19)"""

    q_limit = float(sys.argv[1])
    with open('working/exome.loh_mutations', 'w') as of:
        data_dir = 'data/exome/'
        for afile in global_settings.exome_types:
            file = os.path.join(data_dir, afile)
            write_mutations(file, of, q_limit)

    # with open('working/all_non_ref_hg19.loh_mutations', 'w') as of:
    #     data_dir = 'data/all_non_ref_hg19/'
    #     for exome_type in global_settings.exome_types:
    #         file = os.path.join(data_dir, exome_type)
    #         write_mutations_non_ref(file, of, q_limit)
            
if __name__ == '__main__':
    main()
        
            
