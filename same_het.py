"""Find instance where the same heterozygous call is in cancer and tumor"""
import global_settings, os, sys
from collections import defaultdict
import nose.tools

def init_zero():
    return 0

def write_mutations(afile, of, q_limit):
    """Find instances where the call from normal to cancer changes"""

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
                        if calls[cancer] == calls[normal] and calls[cancer] in global_settings.het_bases:
                            somatics[cancer+':'+normal] += 1
                            of.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'
                                     % (afile.split('/')[-1], chr, pos,
                                        normal, cancer,                                        
                                        calls[normal], calls[cancer]))
                    
    for sample in somatics:
        print afile.split('/')[-1] + '\t' + sample + '\t' + str(float(100)*float(somatics[sample])/float(totals[sample]))

q_limit = float(sys.argv[1])
with open('working/het2het_calls', 'w') as of:
    data_dir = 'data/exome/'
    for afile in os.listdir(data_dir):
        if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
            file = os.path.join(data_dir, afile)
            write_mutations(file, of, q_limit)
            
