"""This holds the parsed data for a sample from
data/all_non_ref_hg19"""
import global_settings, os, utils
from collections import defaultdict

def get_consensus_qualities(afile):
    """Grab quality score for all calls. This is parsing the pileup output (ann files)"""

    qualities = {}
    with open(afile) as f:
        for line in f:
            sp = line.strip().split('\t')
            chr = sp[2]
            pos = sp[3]
            con_quality = float(sp[13])
            qualities[chr+':'+pos] = con_quality
    return qualities

def parse_call_file(call_file, consensus_quality):
    """Grab all info from the file. 
       This only parses the all_non_ref_hg19 files.

       Organization:
       For each chrpos, I have
       ref_allele
       mutation_type (AA:AA)
       snp (in dbSNP if available)
       sample_type (normal/cancer call data, N|T)
            sample_type contains for normal or cancer
            call
            snp_quality
            consensus_quality
            coverage
            call_str
    """

    chrpos_data = {}
    with open(call_file) as f:
        for line in f:
            sp = line.split('\t')
            (gene_name, ref_seq,
             chr, pos, mutation_type,
             snp_name,
             snp, strand,
             sample_count, sample_str,
             something, ref_allele,
             bases_before, bases_after,
             avg_snp_quality,
             min_snp_quality,
             max_snp_quality) = sp[0:17]
            if '_' not in chr and 'M' not in chr:
                chrpos = chr + ':' + pos
                if chrpos not in chrpos_data:
                    chrpos_data[chrpos] = {}
                    chrpos_data[chrpos]['ref_allele'] = ref_allele.upper()
                    chrpos_data[chrpos]['mutation_type'] = mutation_type
                    if snp.strip():
                        chrpos_data[chrpos]['snp'] = snp
                    samples = sample_str.split('-')
                    idx = 18
                    for sample in samples:
                        sample_type = sample[-1] # N|T
                        if sample_type not in ('N', 'T'):
                            raise ValueError
                        chrpos_data[chrpos][sample_type] = {}
                        call = sp[idx].upper()
                        call_snp_quality = float(sp[idx+1])
                        coverage = int(sp[idx+2])
                        call_str = sp[idx+3]
                        for field, val in (('call', call),
                                           ('snp_quality', call_snp_quality),
                                           ('coverage', coverage),
                                           ('call_str', call_str),
                                           ('consensus_quality', 
                                            consensus_quality[sample][chrpos])):
                            chrpos_data[chrpos][sample_type][field] = val
                        idx += 5
    return chrpos_data

def write_MU2A_input(mutations, output_file, file_type):
    """Write MU2A input for a file_type (somatic | inherited)"""

    with open(output_file, 'w') as f:
        for sample in mutations:
            for chrpos in mutations[sample]:
                chr, pos = chrpos.split(':')
                (mutation_type, normal_call, 
                 cancer_call, ref_call) = mutations[sample][chrpos]
                if file_type == 'inherited':
                    f.write(sample + '\t' + chr + '\t' + pos + '\t' 
                            + ref_call + '\t' + cancer_call + '\n')
                elif file_type == 'somatic':
                    f.write(sample + '\t' + chr + '\t' + pos + '\t' 
                            + normal_call + '\t' + cancer_call + '\n')
                else:
                    raise ValueError 
        
class calls:
    """Simple class for data parsed from data/all_non_ref_hg19/yuaker folder, or other samples"""

    def __init__(self, data_dir, sample_name):
        """Parses exome and ann files in data_dir for this sample to make {} of paired tumor/normal call data"""

        self.data = {} # exome_types on top level
        cancer_ann_file = os.path.join(data_dir, sample_name 
                                       + 'T.ann')
        normal_ann_file = os.path.join(data_dir,  sample_name 
                                       + 'N.ann')
        cancer_consensus_qualities = get_consensus_qualities(cancer_ann_file)
        normal_consensus_qualities = get_consensus_qualities(normal_ann_file)
        consensus_qualities = {}
        consensus_qualities[sample_name + 'T'] = cancer_consensus_qualities
        consensus_qualities[sample_name + 'N'] = normal_consensus_qualities
        for exome_type in global_settings.exome_types:
            call_file = os.path.join(data_dir, exome_type)
            self.data[exome_type] = parse_call_file(call_file, consensus_qualities)

    def get_inherited_somatic_mutations(self, quality_cutoff, coverage_cutoff):
        """Return chrpos and ref, normal, cancer alleles for in inherited and somatic {}"""

        # exome types on top level
        inherited = defaultdict(dict)
        somatic = defaultdict(dict)
        for exome_type in self.data:
            for chrpos in self.data[exome_type]:
                if self.data[exome_type][chrpos]['N']['consensus_quality'] > quality_cutoff and self.data[exome_type][chrpos]['T']['consensus_quality'] > quality_cutoff and self.data[exome_type][chrpos]['N']['coverage'] >= coverage_cutoff and  self.data[exome_type][chrpos]['T']['coverage'] >= coverage_cutoff and max((self.data[exome_type][chrpos]['N']['snp_quality'], self.data[exome_type][chrpos]['T']['snp_quality'])) > quality_cutoff:
                    mutation_type = self.data[exome_type][chrpos]['mutation_type']
                    m1, m2 = mutation_type.split(':')
                    if mutation_type in ('BB:BB', 'AB:AB'):
                        inherited[exome_type][chrpos] = (mutation_type,
                                                         self.data[exome_type][chrpos]['N']['call'],
                                                         self.data[exome_type][chrpos]['T']['call'],
                                                         self.data[exome_type][chrpos]['ref_allele'])
                    elif mutation_type not in ('AB:AA', 'AB:BB') and m1 != m2:
                        somatic[exome_type][chrpos] = (mutation_type,
                                                       self.data[exome_type][chrpos]['N']['call'],
                                                       self.data[exome_type][chrpos]['T']['call'],
                                                       self.data[exome_type][chrpos]['ref_allele'])
        # data is organized by exome_type
        return (inherited, somatic)

    def write_mutations(self, quality_cutoff, coverage_cutoff, exome_type, somatic_or_inherited, afile):
        """Write mutations for file"""

        inherited, somatic = self.get_somatic_inherited_mutations(quality_cutoff, coverage_cutoff)
        if somatic_or_inherited == 'somatic':
            use_muts = somatic
        else:
            use_muts = inherited
        with open(afile, 'w') as f:
            for chrpos in inherited[exome_type]:
                f.write(chrpos + '\n')

    def mk_MU2A_input(exome_type, quality_cutoff, coverage_cutoff,
                      inherited_output_file, somatic_output_file):
        """Write input for MU2A for inherited and somatic mutations. This is broken b/c of the function above"""
        
        inherited, somatic = self.get_somatic_inherited_mutations(quality_cutoff, coverage_cutoff)
        write_MU2A_input(inherited[exome_type], inherited_output_file, 'inherited')
        write_MU2A_input(somatic[exome_type], somatic_output_file, 'somatic')

    

def get_mutations_for_paired_samples(quality_cutoff, coverage_cutoff):
    """Load mutations from data/all_non_ref_hg19 for all paired samples"""

    sample2exome2mutations = {}
    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg19/'

    for cancer, normal in global_settings.pairs:
        sample_name = cancer.split('0')[0]
        data_dir = os.path.join(use_data_dir, cancer)
        mutation_calls = calls(data_dir, sample_name)
        (inherited, 
         somatic) = mutation_calls.get_inherited_somatic_mutations(quality_cutoff,
                                                                   coverage_cutoff)
        sample2exome2mutations[sample_name] = (inherited, somatic)

    return sample2exome2mutations

def get_data_for_paired_samples():
    """Grab all data w/ no quality cutoff"""

    sample2data = {}
    use_data_dir = '/home/perry/Projects/loh/data/all_non_ref_hg19/'

    for cancer, normal in global_settings.pairs:
        sample_name = cancer.split('0')[0]
        data_dir = os.path.join(use_data_dir, cancer)
        mutation_calls = calls(data_dir, sample_name)
        sample2data[sample_name] = mutation_calls.data

    return sample2data

def get_coverages(sample2data):
    """For each sample, return chrpos 2 coverage for each exome type.
       Return sample 2 exometype 2 chrpos 2 {N:cov, T:cov}"""

    coverage = {}
    for sample in sample2data:
        coverage[sample] = {}
        for exome_type in sample2data[sample]:
            coverage[sample][exome_type] = {}
            for chrpos in sample2data[sample][exome_type]:
                d = {}
                d['N'] = sample2data[sample][exome_type][chrpos]['N']['coverage']
                d['T'] = sample2data[sample][exome_type][chrpos]['T']['coverage']
                coverage[sample][exome_type][chrpos] = d
    return coverage

def get_normal_val(ref_allele, str_length, str):
    """Count minor allele frequency."""

    ref_count = 0
    base_counts = defaultdict(utils.init_zero)
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

def get_freq_diff(ref_allele, call, str_length, call_str, 
                  normal_freq_and_base):
    """Take diff btwn normal_base_freq and cancer"""

    ref_count = 0
    base_counts = defaultdict(utils.init_zero)
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
    return abs(normal_freq - tumor_freq)

def get_allele_freq_diffs_for_loh(data, quality_cutoff, coverage_cutoff):
    """Find allele frequency differenes for all normal heterozygous"""

    freq_diffs = {}
    for exome_type in data:
        freq_diffs[exome_type] = {}
        for chrpos in data[exome_type]:
            m1, m2 = data[exome_type][chrpos]['mutation_type'].split(':')
             # only take normal heterozygous
            if m1[0] != m1[1] and chrpos not in freq_diffs[exome_type]:

                normal_freq_and_base = get_normal_val(data[exome_type][chrpos]['ref_allele'], 
                                                      data[exome_type][chrpos]['N']['coverage'], 
                                                      data[exome_type][chrpos]['N']['call_str'])
                freq_diffs[exome_type][chrpos] = get_freq_diff(data[exome_type][chrpos]['ref_allele'], 
                                                               data[exome_type][chrpos]['T']['call'], 
                                                               data[exome_type][chrpos]['T']['coverage'], 
                                                               data[exome_type][chrpos]['T']['call_str'], 
                                                               normal_freq_and_base)
    return freq_diffs
        
        
