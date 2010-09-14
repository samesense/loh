"""I'm trying to figure out why we have ~200 more mutations that Murim"""
import os
import cmp_murim_mutations_yusan

def get_my_mutations():
    """Load mutations from working/"""

    my_mutations = {}
    with open('/home/perry/Projects/loh/working/murim.exome.aa_chg.vars') as f:
        for line in f:
            my_mutations[line.strip()] = True
    return my_mutations

def get_murim_mutations():
    """Get coding somatic mutations above 100 quality"""

    os.system("grep 'Coding-missense' /home/perry/Projects/loh/data/murim/CANCER_specific.csv | grep somatic | cut -f 3,4,10 > /home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual")
    murim_mutations = {}
    with open('/home/perry/Projects/loh/data/murim/somatic_coding_parsed_chr_pos_qual') as f:
        for line in f:
            sp = line.strip().split('\t')
            if float(sp[-1]) > float(100): # murim uses 80
                murim_mutations[sp[0] + ':' + sp[1]] = True
    return murim_mutations

#print 'both', len(set(my_mutations) & set(murim_mutations))
#print 'us', len(my_mutations)
#print 'murim', len(murim_mutations)

def get_mutations_we_miss():
    """Return chrpos not in our results"""

    murim_mutations = get_murim_mutations()
    my_mutations = get_my_mutations()

    missing = {}
    for chrpos in murim_mutations:
        if chrpos not in my_mutations:
            missing[chrpos] = True
    return missing
    
