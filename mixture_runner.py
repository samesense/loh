"""Make CNV/LOH plots and mixture plot"""
import call_class, os, random, global_settings, utils

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
    sample_pairs = []
    for sample in coverages:
        sample_pairs.append((sample + 'T',
                             sample + 'N'))
    return sample_pairs

def CNV_plot(cnv_file, subdir, png):
    """Make CNV plots using R cnv library"""

    os.system('mkdir -p working/cnv_seq/CNV')
    os.system('mkdir -p working/cnv_seq/CNV/' 
              + subdir + '/')
    output_cnv_file = 'working/cnv_seq/CNV/' + subdir + '/' + cnv_file.split('.hits')[0]  + '.cnvs'
    os.system('rm ' + output_cnv_file)
    rtmp = 'rtmp' + str(random.randint(0,1000))
    with open(rtmp, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write('library(cnv)\n')
        f.write("data<-read.delim('"
                + cnv_file + "')\n")
        f.write("png('" + png + "')\n")
        f.write("plot.cnv.all.perry(data,colour=9)\n")
        f.write('dev.off()\n')
        f.write("cnv.print(data, file='" + output_cnv_file + "')\n")
        f.write('q()\n')
    if utils.check_input(cnv_file):
        os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')
        os.system('rm ' + rtmp + ' tmpLog')
    os.system('mv ' + cnv_file + ' working/cnv_seq/CNV/' + subdir + '/')
    os.system('mv ' + cnv_file.replace('cnv', 'count') + ' working/cnv_seq/CNV/' + subdir + '/')

def call_cnv_seq(sample_pairs):
    """System calls to CNV-seq"""

    return_plots = {}
    subdir = 'all_non_ref_hg19'
    os.system('mkdir -p ' + os.path.join('plots', 'cnv-seq'))
    os.system('mkdir -p ' + os.path.join('plots', 'cnv-seq', subdir))
    for exome in global_settings.exome_types:
        for cancer, normal in sample_pairs:
            cancer_hits = 'working/cnv_seq/' + subdir + '/' + exome + '.'  + cancer + '.hits'
            normal_hits = 'working/cnv_seq/' + subdir + '/' + exome + '.' + normal + '.hits'
            os.system('perl cnv-seq.pl --test ' + cancer_hits + ' --ref '
                      + normal_hits 
                      + ' --genome human --log2 0.6 --p 0.001 --bigger-window 1.5 --annotate -minimum-windows 4')
            cnv_file = exome + '.' + cancer + '.hits-vs-' + exome + '.' + normal + '.hits.log2-0.6.pvalue-0.001.minw-4.cnv'
            plot_file = 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer + '.png'
            CNV_plot(cnv_file, subdir, plot_file)  
            return_plots[cancer[0:-1]] = plot_file
    return return_plots

def plot_loh():
    """Recreate Murim's plots"""

    os.system('mkdir -p working/loh/')
    samples2data = call_class.get_data_for_paired_samples()
    

def main():
    """Entry point"""

    random.seed()
    sample_pairs = mk_cnv_seq_input_runner()
    cnv_plots = call_cnv_seq(sample_pairs)
    loh_plots = plot_loh()
    for sample in cnv_plots:
         os.system('montage -geometry 500 -quality 100 -tile 1x2 '
                   cnv_plots[sample] + ' '
                   loh_plots[sample] + ' '
                   + cnv_plots[sample].replace('png', 'cnv_loh.png')

if __name__ == '__main__':
    main()
