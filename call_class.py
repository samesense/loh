"""This holds the parsed data for a sample from
   data/all_non_ref_hg19"""
import global_settings, os
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

    def get_somatic_inherited_mutations(self, quality_cutoff, coverage_cutoff):
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
