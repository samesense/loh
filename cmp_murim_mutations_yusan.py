"""Compare our yusan normal->cancer mutations
   to Murim's yusan normal->cancer mutations.
   This work is done w/ hg18."""
import os

def get_mutations(normal, cancer):
    """Given calls in normal and cancer, return the differences"""

    # chrpos
    mutants = {}
    seen_chrpos = {}
    for chrpos in normal:
        if chrpos in cancer:
            if normal[chrpos] != cancer[chrpos]:
                mutants[chrpos] = normal[chrpos] + ':' + cancer[chrpos]
            seen_chrpos[chrpos] = True
    return (mutants, seen_chrpos)

def get_mutations_limited(normal, cancer, limited):
    """Given calls in normal and cancer, return the differences"""

    # chrpos
    mutants = {}
    for chrpos in normal:
        if chrpos in cancer and chrpos in limited:
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
            have_data = False
            if line[0] != '#':
                sp = line.split('\t')
                
                # snpfilter
#                 call = sp[5]
#                 quality = sp[9]
#                 chr = sp[2]
#                 base = sp[3]
#                have_data = True
                
                # shortall
                if line[0:3] == 'chr':
                    call = sp[3]
                    quality = sp[6]
                    chr = sp[0]
                    base = sp[1]
                    have_data = True
                if have_data:
                    if int(quality) > quality_cutoff:
                        calls[chr+':'+base] = call
    return calls

quality_cut = 100
murim_data_dir = 'data/murim/'
#murim_suffix_normal = 'EX_BLD1_1ln.snpfilter_anno_v1.62.txt'
#murim_suffix_cancer = 'EX_CANC_1ln.snpfilter_anno_v1.62.txt'
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
us_mutants, us_chrpos = get_mutations(us_normal, us_cancer)

print 'US normal/cancer mutants'
print str(len(us_normal)) + '/' + str(len(us_cancer)), len(us_mutants)

murim_normal = load_murims_calls(os.path.join(murim_data_dir,
                                              'yusanN',
                                              murim_suffix_normal), 
                                 quality_cut)
murim_cancer = load_murims_calls(os.path.join(murim_data_dir,
                                              'yusanT',
                                              murim_suffix_cancer), 
                                 quality_cut)
murim_mutants, murim_chrpos = get_mutations(murim_normal, murim_cancer)

murim_limited_mutants = get_mutations_limited(murim_normal, murim_cancer, us_chrpos)
us_limited_mutants = get_mutations_limited(us_normal, us_cancer, murim_chrpos)

print 'murim normal/cancer mutants'
print str(len(murim_normal)) + '/' +  str(len(murim_cancer)), len(murim_mutants)

print 'us normal, murim normal differences', len(get_mutations(murim_normal, us_normal))
print 'final mutation estimates us/murim', str(len(us_limited_mutants))+ '/'+ str(len(murim_limited_mutants))


