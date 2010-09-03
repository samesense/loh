"""Compare our yusan normal->cancer mutations
   to Murim's yusan normal->cancer mutations.
   This work is done w/ hg18."""
import os

def get_mutations(normal, cancer):
    """Given calls in normal and cancer, return the differences"""

    # chrpos
    mutants = {}
    for chrpos in normal:
        if chrpos in cancer:
            if normal[chrpos] != cancer[chrpos]:
                mutants[chrpos] = normal[chrpos] + ':' + cancer[chrpos]
    return mutants

def load_our_calls(afile, quality_cutoff):
    """Load our differenes from ref w/ some quality 
       for the 1st # (read quality"""

    # chrpos: call
    calls = {}
    with open(afile) as f:
        for line in f:
            if line.split('\t')[2] != '*':
                (chr, base, ref, call, read_qual, 
                 notsure1, notsure2, reads, read_str,
                 something) = line.strip().split('\t')
                if int(read_qual) > quality_cutoff and int(reads) >= 8:
                    calls[chr + ':' + base] = call
    return calls

def load_murims_calls(afile, quality_cutoff):
    """Load Murim's differences from ref w/ some quality.
       This parses the shortall file."""

    # chrpos: call
    calls = {}
    with open(afile) as f:
        for line in f:
            sp = line.split('\t')
            if len(sp) == 29:
                if '#' not in sp[0]:
                    call = sp[3]
                    quality = sp[6]
                    chr = sp[0]
                    base = sp[1]
                    calls[chr+':'+base] = base
    return calls

quality_cut = 100
murim_data_dir = 'data/murim/'
murim_suffix_normal = 'EX_BLD1_1ln.snpfilter_anno_shortall_v1.62.txt'
murim_suffix_cancer = 'EX_CANC_1ln.snpfilter_anno_shortall_v1.62.txt'
our_data_dir = 'data/raw_hg18/'
our_suffix = 's_5.var'

us_normal = load_our_calls(os.path.join(our_data_dir,
                                        'yusanN', our_suffix),
                           quality_cut)
us_cancer = load_our_calls(os.path.join(our_data_dir,
                                        'yusanT', our_suffix),
                           quality_cut)
us_mutants = get_mutations(us_normal, us_cancer)
print 'US normal/cancer mutants'
print len(us_normal), len(us_cancer), len(us_mutants)

murim_normal = load_murims_calls(os.path.join(murim_data_dir,
                                              'yusanN',
                                              murim_suffix_normal), 
                                 quality_cut)
murim_cancer = load_murims_calls(os.path.join(murim_data_dir,
                                              'yusanT',
                                              murim_suffix_cancer), 
                                 quality_cut)
murim_mutants = get_mutations(murim_normal, murim_cancer)

print 'murim normal/cancer mutants'
print len(murim_normal), len(murim_cancer), len(murim_mutants)


