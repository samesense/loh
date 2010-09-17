"""Make CNV/LOH plots and mixture plot"""
import call_class, os

def mk_cnv_seq_input_chrpos(coverage, chr, pos, file):
    """ Print out each position the # of times it is covered (signal strength)."""

    for i in xrange(coverage):
        file.write('%s\t%s\n' % (chr.split('chr')[1], pos))

def mk_cnv_seq_input(coverages, output_dir):
    """Make input for CNV-seq.
       Print out each position the # of times it is covered (signal strength)."""

    for sample in coverages:
        for exome_type in coverages[sample]:
            with open(os.path.join(output_dir,
                                   '.'.join([exome_type,
                                            sample + 'N.hits'])), 'w') as nhits:
                with open(os.path.join(output_dir,
                                       '.'.join([exome_type,
                                                sample + 'T.hits'])), 'w') as thits:
                    for chrpos in coverages[sample][exome_type]:
                        chr, pos = chrpos.split(':')
                        for coverage, file in ((coverages[sample][exome_type][chrpos]['N'],
                                                nhits),
                                               (coverages[sample][exome_type][chrpos]['T'],
                                                thits)):
                            mk_cnv_seq_input_chrpos(coverage, chr, pos, file)

def mk_cnv_seq_input_runner():
    """Construct CNV-seq input from all_non_ref_hg19 5 paired samples"""

    full_data_dir = 'all_non_ref_hg19/'
    cnv_seq_dir = 'working/cnv_seq/'
    os.system('mkdir -p ' + cnv_seq_dir)
    output_dir = os.path.join(cnv_seq_dir, full_data_dir)
    os.system('mkdir -p ' + output_dir)
    samples2data = call_class.get_data_for_paired_samples()
    coverages = call_class.get_coverages(samples2data)
    mk_cnv_seq_input(coverages, output_dir)

def main():
    """Entry point"""

    mk_cnv_seq_input_runner()
    
if __name__ == '__main__':
    main()
