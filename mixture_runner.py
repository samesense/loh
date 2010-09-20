"""Make CNV/LOH plots and mixture plot"""
import call_class, os, random, global_settings, utils, sys
from collections import defaultdict

def avg_diffs_per_chr(allele_diffs, cnvs, mix_file):
    """For each chr, use the lost of copy CNVs to compuete the averag allele diff per chr"""
    
    with open(mix_file, 'w') as f:
        for chr in allele_diffs:
            sum_diffs = float(0)
            total_calls = 0
            c = 'chr' + chr
            if c in cnvs:
                for pos in allele_diffs[chr]:
                    for cnv_st, cnv_end in cnvs[c]:
                        if pos < cnv_st:
                            break
                        if pos >= cnv_st and pos <= cnv_end:
                            sum_diffs += allele_diffs[chr][pos]
                            total_calls += 1
                            break
            if total_calls:
                avg_diff = str(sum_diffs/float(total_calls))
                f.write('%s\t%d\t%s\n' %
                        (chr, total_calls, avg_diff))
            
def get_cnvs(afile):
    """Load los of copy CNVs computed by cnv-seq"""

    cnvs = defaultdict(list)
    with open(afile) as f:
        f.readline()
        for line in f:
            sp = line.strip().split('\t')
            if len(sp) == 7:
                cnv, chr, st, end, size, log2, pval = line.strip().split('\t')
                if log2 != 'NA':
                    # check for copy loss
                    if float(log2) < float(-.8):
                        cnvs[chr].append((int(st), int(end)))
    for chr in cnvs:
        cnvs[chr].sort()
    return cnvs

def get_allele_diffs(afile):
    """Load Murim's normal/cancer diffs for heterozygous calls in normal."""

    allele_diffs = defaultdict(dict)
    with open(afile) as f:
        f.readline()
        for line in f:
            chr, pos, diff = line.strip().split('\t')
            allele_diffs[chr][int(pos)] = float(diff)
    return allele_diffs

def find_mixture(loh_diff_file, cnv_file, mix_file):
    """This is the cumulation of the CNV/LOH work. We need the mixture of normal/cancer cells to help call somatic mutations. I will find the mixture per chromosome by looking at copy loss regions, as indicated by CNV-seq, and then look at the allele differences (Murim's plot) to determine the mixture by looking at the difference between 0.5 and the observed normal/cancer allele frequency differences."""

    allele_diffs = get_allele_diffs(loh_diff_file)
    cnvs = get_cnvs(cnv_file)
    avg_diffs_per_chr(allele_diffs, cnvs, mix_file)

def plot_mixture(snp_cut):
    """Use R ggplot2 to plot the mixture for all exome types and cancer/normal pairs"""

    subdir = 'all_non_ref_hg19'
    os.system('mkdir -p working/mixture/')
    os.system('mkdir -p working/mixture/all_non_ref_hg19/')
    for exome_type in global_settings.exome_types:
        for cancer, normal in global_settings.pairs:
            cnv_file = 'working/cnv_seq/CNV/' + subdir + '/' + exome_type + '.' + cancer.split('0')[0].strip('T') + '.cnvs'
            loh_file = 'working/loh/' + subdir + '/' + exome_type + '.' + cancer.split('0')[0].strip('T')
            mix_file = 'working/mixture/' + subdir + '/' + exome_type + '.' + cancer.split('0')[0].strip('T') + '.mix'
            find_mixture(loh_file, cnv_file, mix_file)
    rinput = 'rinput' + str(random.randint(0,1000))
    working_dir = os.path.join('working/mixture/', subdir + '/')
    with open(rinput, 'w') as f:
        f.write('Sample\tExome\tChr\tMixture\n')
        for exome_type in global_settings.exome_types:
            for cancer, normal in global_settings.pairs:
                mix_file_name = working_dir + exome_type + '.' + cancer + '.mix'
                if os.path.exists(mix_file_name):
                    with open(mix_file_name) as mixfile:
                        for line in mixfile:
                            chr, snp_count, avg_diff = line.strip().split('\t')
                            if int(snp_count) > snp_cut:
                                f.write('%s\t%s\t%s\t%f\n' %
                                        (cancer, exome_type.split('.')[1], chr, 
                                         float(0.5)-float(avg_diff)))

    # make R calls
    rtmp = 'rtmp' + str(random.randint(0,1000))
    with open(rtmp, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("data<-read.delim('" + rinput + "',header=TRUE,sep='\\t')\n")
        f.write("png('plots/" + subdir + ".mixture.png')\n")
        f.write('ggplot(data) + aes(x=Chr,y=Mixture) + geom_point() + facet_grid(Sample~Exome)\n')
        f.write('dev.off()\n')
        f.write('q()\n')

    os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')
    os.system('rm ' + rinput + ' ' + rtmp)

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
    output_cnv_file = 'working/cnv_seq/CNV/' + subdir + '/' + cnv_file.split('.hits')[0].rstrip('T')  + '.cnvs'
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
    os.system('mv ' + cnv_file.replace('cnv', 'count') 
              + ' working/cnv_seq/CNV/' + subdir + '/')

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
            plot_file = 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer.rstrip('T') + '.png'
            CNV_plot(cnv_file, subdir, plot_file)  
            return_plots[cancer[0:-1] + '.' + exome] = plot_file
    return return_plots

def plot_loh(freq_diffs, plot_file, input_file):
    """Recreate Murim's plots"""

    with open(input_file, 'w') as f:
        f.write('CHR\tMapInfo\tDiff\n')    
        for chr_pos in freq_diffs:
            chr, pos = chr_pos.split(':')
            f.write(chr.split('chr')[1] + '\t' + pos + '\t' + str(freq_diffs[chr_pos]) + '\n')
            
    rtmp = 'rtmp' + str(random.randint(0,1000))
    with open(rtmp, 'w') as f:
        f.write("source('funcs.R')\n")
        f.write("data<-read.delim('"
                + input_file + "')\n")
        f.write("png('" + plot_file + "')\n")
        f.write("plot.murim(data,colour=9)\n")
        f.write('dev.off()\n')
        f.write('q()\n')

    os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')
    os.system('rm ' + rtmp + ' tmpLog')

def loh():
    """Make data for Murim's plots and call function to plot them"""

    plots = {}
    os.system('mkdir -p working/loh/')
    os.system('mkdir -p working/loh/all_non_ref_hg19')
    os.system('mkdir -p plots/loh/')
    os.system('mkdir -p plots/loh/all_non_ref_hg19/')
    quality_cutoff = float(100)
    coverage_cutoff = 8
    samples2data = call_class.get_data_for_paired_samples()
    for sample in samples2data:
        freq_diffs = call_class.get_allele_freq_diffs_for_loh(samples2data[sample], 
                                                              quality_cutoff, 
                                                              coverage_cutoff)
        for exome in freq_diffs:
            plot_file = 'plots/loh/all_non_ref_hg19/' + exome + '.' + sample + '.png'
            input_file = 'working/loh/all_non_ref_hg19/' + exome + '.' + sample
            plot_loh(freq_diffs[exome], plot_file, input_file)
            plots[sample + '.' + exome] = plot_file
    return plots

def main():
    """Entry point"""

    snp_cut = int(sys.argv[1])
    random.seed()
    sample_pairs = mk_cnv_seq_input_runner()
    cnv_plots = call_cnv_seq(sample_pairs)
    loh_plots = loh()
    for sample_exome in cnv_plots:
         os.system('montage -geometry 500 -quality 100 -tile 1x2 '
                   + cnv_plots[sample_exome] + ' '
                   + loh_plots[sample_exome] + ' '
                   + cnv_plots[sample_exome].replace('png', 'cnv_loh.png'))
    plot_mixture(snp_cut)

if __name__ == '__main__':
    main()
