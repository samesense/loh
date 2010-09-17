"""This holds the parsed data for a sample from
   data/all_non_ref_hg19"""
import global_settings

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

def parse_call_file(call_file, data_to_fill):
    """Grab all info from the file.
       Fields """

class calls:
    """Simple class for data parsed from data/all_non_ref_hg19/yuaker folder, or other samples"""

    def __init__(self, data_dir, quality_cutoff, coverage_cutoff):
        """Parses exome and ann files in data_dir for this sample"""

        self.data = {}
        sample_name = data_dir.split('/')[-1]
        cancer_ann_file = os.path.join(data_dir, sample_name 
                                       + 'T.ann')
        normal_ann_file = os.path.join(data_dir,  sample_name 
                                       + 'N.ann')
        cancer_qualities = get_consensus_qualities(cancer_ann_file)
        normal_qualities = get_consensus_qualities(normal_ann_file)
        for exome in global_settings.exome_types:
            call_file = os.path.join(data_dir, exome_type)
            parse_call_file(call_file, self.data,
                            quality_cutoff, coverage_cutoff,
                            cancer_qualities, normal_qualities)
