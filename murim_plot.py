"""Recreate Murim's plot for paired normal/cancer samples.
   Choose to plot all normal heterozygous OR
   Only use normal->cancer het->homo calls."""
import os, sys, global_settings, mutations
from collections import defaultdict

# cancer, normal
pairs = global_settings.pairs
cancer2normal = {}
for cancer, normal in pairs:
    cancer2normal[cancer] = normal
normals = [x[1] for x in pairs]
cancers = [x[0] for x in pairs]
all = []
for n in normals: all.append(n)
for c in cancers: all.append(c)

def get_limiting_mutations(limiting_file):
    """Load normal->cancer het->homo mutations found by loh.py.
       Make a {} of exome_type to chr:pos to normal samples.
       Also load het->het (R,R) | (Y,Y) etc.
       Limit the mutations I'm examining in some way."""

    exome2mutations = {}
    with open(limiting_file) as f:
        for line in f:
            exome, chr, pos, normal, cancer, n_call, c_call = line.strip().split('\t')
            if not exome in exome2mutations:
                exome2mutations[exome] = {}
            if not normal in exome2mutations[exome]:
                exome2mutations[exome][normal] = {}
            loc = chr+':'+pos
            exome2mutations[exome][normal][loc] = True

    return exome2mutations

def init_zero():
    return 0

def get_normal_val(ref_allele, str_length, str):
    """Count minor allele frequency."""

    ref_count = 0
    base_counts = defaultdict(init_zero)
    for s in str:
        if s == '.' or s == ',':
            ref_count += 1
        elif s.upper() in global_settings.homo_bases:
            base_counts[s.upper()] += 1
    base_count_ls = []
    for base in base_counts:
        base_count_ls.append((base_counts[base], base))
    base_count_ls.sort()
    if ref_allele == base_count_ls[-1][1]:
        use_base = base_count_ls[-2][1]
    else:
        use_base = base_count_ls[-1][1]
        
    return (float(base_counts[use_base])/float(str_length), use_base)

def get_het_normal(file):
    """Grab hetero SNPs from normal samples"""
    
    normal2het = defaultdict(dict)
    with open(file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            for s in samples:
                if s in normals:
                    ref_allele = sp[10]
                    chr = sp[2].split('chr')[1]
                    pos = sp[3]
                    quality = float(sp[idx+1])
                    coverage = int(sp[idx+2])
                    call = sp[idx]
                    # use quality > 100
                    # keep if using at least 8x coverage
                    # and if the chr is standard
                    # i.e. no 6_cox_hap2
                    if call not in global_settings.homo_bases and coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
                       #  if 'XY' == chr:
#                             chr = '100'
#                         elif 'X' == chr:
#                             chr = '98'
#                         elif 'Y' == chr:
#                             chr = '99' 
                        normal2het[s][chr+':'+pos] = get_normal_val(ref_allele,
                                                                    coverage,
                                                                    sp[idx+3])
                idx += 5
    return normal2het

def get_norm_sample(sample):
    """If this sample is cancer, and pairs
       with a normal, return that normal.
       Otherwise return False."""

    if sample in cancer2normal:
        return cancer2normal[sample]
    else:
        return False

def get_freq(ref_allele, call, str_length, call_str, normal_freq_and_base, chr, pos, normal_name, exome_type):
    """Take diff btwn normal_base_freq and cancer"""

    ref_count = 0
    base_counts = defaultdict(init_zero)
    for s in call_str:
        if s == '.' or s == ',':
            ref_count += 1
        elif s.upper() in global_settings.homo_bases:
            base_counts[s.upper()] += 1
    normal_freq, normal_base = normal_freq_and_base
    
    tumor_freq = float(0)
    if normal_base in base_counts:
        tumor_freq = float(base_counts[normal_base])/float(str_length)
    elif normal_base == call:
        tumor_freq = float(ref_count)/float(str_length)
#    print str(normal_freq) + '\t' + str(tumor_freq) + '\t' + chr + '\t' + pos + '\t' + normal_name + '\t' + exome_type
    return abs(normal_freq - tumor_freq)

def print_control(chr_pos, normal_freq_base, o):
    """Write out data for control(blue/normal) plots"""

    chr, pos = chr_pos.split(':')
    freq, base = normal_freq_base
    allele_freq = abs(float(.5) - freq)
    o.write(chr + '\t' + pos + '\t' 
            + str(allele_freq) + '\n')

def eval_cancer_diff(sp, chr, pos, normal2het, norm_sample, sample2alleles, tumor_sample, idx, afile):
    """Record cancer/normal difference"""

    quality = float(sp[idx+1])
    coverage = int(sp[idx+2])
    ref_allele = sp[10]
    # use quality > 100
    # keep if using at least 8x coverage
    # and if the chr is standard
    # i.e. no 6_cox_hap2
    # sample must be het in normal or paired w/ het in normal
    if chr+':'+pos in normal2het[norm_sample] and coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
        # if 'XY' == chr:
#             chr = '100'
#         elif 'X' == chr:
#             chr = '98'
#         elif 'Y' == chr:
#             chr = '99'
        if chr + ':' + pos not in sample2alleles[tumor_sample]:
            sample2alleles[tumor_sample][chr + ':' + pos] = get_freq(ref_allele,
                                                                     sp[idx],
                                                                     coverage,
                                                                     sp[idx+3],
                                                                     normal2het[norm_sample][chr+':'+pos],
                                                                     chr, pos, norm_sample,
                                                                     afile.split('/')[-1])

def do_exome(use_homo, limiting_mutation_file):
    """Write out exome mutations"""

    # this method relies on cancer call being
    # different from the ref
    # I need to add cases where it is the same
    # and thus not present in the file

    os.system('mkdir -p working/murim_plot/exome/')
    
    if use_homo == 'homo':
        exome2mutations = get_limiting_mutations(limiting_mutation_file)
    data_dir = 'data/exome/'
    for afile in global_settings.exome_types:
        file = os.path.join(data_dir, afile)
        normal2het = get_het_normal(file)
        sample2alleles = defaultdict(dict)
        with open(file) as f:
            for line in f:
                sp = line.split('\t')
                samples = sp[8].split('-')
                idx = 17
                chr = sp[2].split('chr')[1]
                pos = sp[3]
                for s in samples:
                    norm_sample = get_norm_sample(s)
                    if norm_sample:
                        if use_homo == 'homo':
                            if chr+':'+pos in exome2mutations[afile][norm_sample]:
                                eval_cancer_diff(sp, chr, pos, normal2het, norm_sample, 
                                                 sample2alleles, s, idx, afile)
                        else:
                            eval_cancer_diff(sp, chr, pos, normal2het, norm_sample, 
                                             sample2alleles, s, idx, afile)
                    idx += 5
        outdir = os.path.join('working', 'murim_plot', 'exome',
                              afile.replace('.', '_'))
        os.system('mkdir -p ' + outdir)
        for s in sample2alleles:
            ofile = os.path.join(outdir, 'hg19_murim.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tDiff\n')
                for chr_pos in sample2alleles[s]:
                    chr, pos = chr_pos.split(':')
                    allele_freq = sample2alleles[s][chr_pos]
                    o.write(chr + '\t' + pos + '\t' + str(allele_freq) + '\n')
        for s in normal2het:
            ofile = os.path.join(outdir, 'hg19_murim.' + s)
            with open(ofile, 'w') as o:
                o.write('CHR\tMapInfo\tDiff\n')
                for chr_pos in normal2het[s]:
                    if use_homo == 'homo':
                        if chr_pos in exome2mutations[afile][s]:
                            print_control(chr_pos, normal2het[s][chr_pos], o)
                    else:
                        print_control(chr_pos, normal2het[s][chr_pos], o)

def do_all_non_ref_freqs(data_file, normal_qualities, cancer_qualities, quality_cutoff, coverage_cutoff):
    """Find allele frequency difference for all_non_ref data"""

    chr2freqDiff = {}
    with open(data_file) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[9].split('-')
            ref_allele = sp[11]
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
            normal_call_str = sp[21]
            cancer_call = sp[23].upper()
            cancer_call_str = sp[26]
            
            # throw exception if normal/cancer name does not match
            if sample_name_normal != sample_name_cancer:
                i = 1/0

            (avg_snp_quality, min_snp_quality, 
             max_snp_quality) = (float(x) for x in sp[14:17])
            if normal_qualities[chrpos] > quality_cutoff and cancer_qualities[chrpos] > quality_cutoff and normal_coverage >= coverage_cutoff and cancer_coverage >= coverage_cutoff and max_snp_quality > quality_cutoff and '_' not in chr and 'M' not in chr:
                m1, m2 = mutation_type.split(':')
                # only take normal heterozygous
                if m1[0] != m1[1] and chrpos not in chr2freqDiff:
                    # call for normal
                    normal_freq_and_base = get_normal_val(ref_allele, normal_coverage, normal_call_str)
                    # call for cancer
                    freq_diff = get_freq(ref_allele,
                                         cancer_call,
                                         cancer_coverage,
                                         cancer_call_str,
                                         normal_freq_and_base,
                                         chr, pos, normal,
                                         data_file.split('/')[-1])
                    chr2freqDiff[chrpos] = freq_diff    
    return chr2freqDiff

def do_all_non_ref():
    """This works with the full dataset and captures more mutations"""

    quality_cutoff = float(100)
    coverage_cutoff = 8
    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg19/'
    os.system('mkdir -p working/murim_plot/all_non_ref_hg19/')
    for exome_type in global_settings.exome_types:
        outdir = os.path.join('working', 'murim_plot', 'all_non_ref_hg19',
                              exome_type.replace('.', '_'))
        os.system('mkdir -p ' + outdir)
    for cancer, normal in global_settings.pairs:
        data_subdir = os.path.join(use_data_dir, cancer)
        sample_name = cancer.split('0')[0]
        cancer_qualities = mutations.get_consensus_qualities(os.path.join(data_subdir,
                                                                          sample_name + 'T.ann'))
        normal_qualities = mutations.get_consensus_qualities(os.path.join(data_subdir,
                                                                          sample_name + 'N.ann'))
        for exome in global_settings.exome_types:
            data_file = os.path.join(data_subdir, exome)
            chrpos2diff = do_all_non_ref_freqs(data_file, 
                                               normal_qualities, 
                                               cancer_qualities, 
                                               quality_cutoff, 
                                               coverage_cutoff)
            with open(os.path.join('working/murim_plot/all_non_ref_hg19/',
                                   exome.replace('.', '_'),
                                   'hg19_murim.' + cancer), 'w') as outf:
                outf.write('CHR\tMapInfo\tDiff\n')
                for chrpos in chrpos2diff:
                    # print 'here', os.path.join('working/murim_plot/all_non_ref_hg19/',
                    #                exome_type.replace('.', '_'),
                    #                'hg19_murim.' + cancer)
                    chr, pos = chrpos.split(':')
                    outf.write(chr.split('chr')[1] + '\t' + pos + '\t' + str(chrpos2diff[chrpos]) + '\n')
    
def main():  
    """Script entry"""

    use_homo = sys.argv[1] # homo | all
    limiting_mutation_file = sys.argv[2]
    # if use_homo == all, limiting_mutation_file does not matter
    os.system('mkdir -p working/murim_plot')
    #do_exome(use_homo, limiting_mutation_file) # this is for verification
    do_all_non_ref() # this is the final product

if __name__ == '__main__':
    main()


