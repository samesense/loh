"""Find somatic changes in cancer."""
import global_settings, os
from collections import defaultdict

def init_zero():
    return 0

def is_loh(cancer, normal):
    """Return true if normal->cancer is LOH"""
    
    return False

def write_mutations(afile, of):
    """Find instances where the call from normal to cancer changes"""

    somatics = defaultdict(init_zero)
    totals = defaultdict(init_zero)
    with open(file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            calls = {}
            for s in samples:
                chr = sp[2].split('chr')[1]
                pos = sp[3]
                quality = float(sp[idx+1])
                coverage = int(sp[idx+2])
                call = sp[idx]
                calls[s] = call
                idx += 5
            # check for changes in normal/cancer pairs
            for cancer, normal in global_settings.pairs:
                if cancer in calls and normal in calls:
                    totals[cancer+':'+normal] += 1
                    if calls[cancer] != calls[normal]:
                        if not is_loh(calls[cancer], calls[normal]):
                            somatics[cancer+':'+normal] += 1
                            of.write('%s\t%s\t%s\t%s\t%s\n'
                                     % (afile.split('/')[-1], normal, cancer, 
                                        calls[normal], calls[cancer]))
                    
    for sample in somatics:
        print 'total', afile.split('/')[-1], sample, float(100)*float(somatics[sample])/float(totals[sample])

with open('working/somatic_mutations', 'w') as of:
    data_dir = 'data/exome/'
    for afile in os.listdir(data_dir):
        if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
            file = os.path.join(data_dir, afile)
            write_mutations(file, of)
            
