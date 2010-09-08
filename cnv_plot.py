"""Make the copy number plot from
   High-resolution genomic profiling
   of chromosomal abberations using
   Infinium whole-genome genotyping.
   Make a plot for normal and cancer.

   Plot the ratio of coverage (not freq)
   for normal and cancer.
"""
import os, random, sys, math, global_settings
from collections import defaultdict

random.seed()

def mk_cvn_seq(coverages, cancer, normal, cancer_file, normal_file):
    """Make input to CNV-seq. Print out each position the # of times it is covered (signal strength)."""

    with open(cancer_file, 'w') as cancerf:
        with open(normal_file, 'w') as normalf:
            for chrpos in coverages[cancer]:
                chr, pos = chrpos.split(':')
                if chrpos in coverages[normal]:
                    chr = chr.split('chr')[1]
                    for i in xrange(coverages[cancer][chrpos]):
                        cancerf.write('%s\t%s\n' %
                                      (chr, pos))
                    for i in xrange(coverages[normal][chrpos]):
                        normalf.write('%s\t%s\n' %
                                      (chr, pos))

def mk_signal_ratio(coverages, cancer, normal, input_file):
    """Write the ratio of cancer:normal signal strength for positions."""

    with open(input_file, 'w') as f:
        f.write('Chr\tPos\tCovRatio\n')
        for chrpos in coverages[cancer]:
            chr, pos = chrpos.split(':')
            if chrpos in coverages[normal]:
                if chr in ['chr1', 'chr6', 'chr2', 'chr10']:
                    f.write('%s\t%s\t%f\n' %
                            (chr, pos, math.log(float(coverages[cancer][chrpos])/
                                                float(coverages[normal][chrpos]), 2)))

def mk_r_signal_ratio_file(r_file, input_file, plot_file):
    """Make R input commands for signal ratio"""

    with open(r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("data<-read.delim('" + input_file 
                + "',header=TRUE,sep='\\t')\n")
        f.write("png('" + plot_file + "')\n")
        f.write("ggplot(data) + aes(x=Pos,y=CovRatio) + geom_point() + facet_grid(Chr~.) + opts(legend.position='none',title='Log Coverage Ratio')\n")
        f.write('dev.off()\n')
        f.write('q()\n')
        
    os.system('R --vanilla < ' + r_file)

def plot_coverage_ratio(coverages, cancer, normal, plot_file):
    """For each chr, plot the coverage ratio of cancer to normal along the chr"""

    input_file = 'input' + str(random.randint(0,1000))
    rtmp = 'rtmp' + str(random.randint(0,1000))

    mk_signal_ratio(coverages, cancer, normal, input_file)
    mk_r_signal_ratio_file(rtmp, input_file, plot_file)
    os.system('rm ' + input_file + ' ' + rtmp)

def mk_r_ratio_input(chrpos2alleleFreqs, input_file):
    """Write our chr pos and ratio of ref allele freq to other allele freq input_file"""

    with open(input_file, 'w') as f:
        f.write('Chr\tPos\tLogRatio\n')
        for chrpos in chrpos2alleleFreqs:
            chr, pos = chrpos.split(':')
            ref_allele_freq, min_allele_freq = chrpos2alleleFreqs[chrpos]
            #if chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
            if ref_allele_freq != float(0):
                f.write('%s\t%s\t%f\n' %
                        (chr, pos, math.log(ref_allele_freq/min_allele_freq)))

def mk_r_input(chrpos2alleleFreqs, input_file):
    """Write our chr pos and allele freq to input_file"""

    with open(input_file, 'w') as f:
        f.write('Chr\tPos\tRefFreq\n')
        for chrpos in chrpos2alleleFreqs:
            chr, pos = chrpos.split(':')
            ref_allele_freq, min_allele_freq = chrpos2alleleFreqs[chrpos]
            #if chr in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']:
            f.write('%s\t%s\t%f\n' %
                    (chr, pos, ref_allele_freq))

def mk_r_file(input_file, r_file, plot_file):
    """Write out R commands for allele freq plot"""
    
    # opts(legend.position='none', axis.text.y = theme_blank(), axis.text.x = theme_blank()
    with open(r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("data<-read.delim('" + input_file 
                + "',header=TRUE,sep='\\t')\n")
        f.write("png('" + plot_file + "')\n")
        f.write("ggplot(data) + aes(x=Pos,y=RefFreq) + geom_point() + facet_grid(Chr~.) + opts(legend.position='none',title='Ref Freq') + ylim(c(0,1))\n")
        f.write('dev.off()\n')
        
    os.system('R --vanilla < ' + r_file)

def mk_r_ratio_file(input_file, r_file, plot_file):
    """Write out R commands for allele ratio plot"""
    
    # opts(legend.position='none', axis.text.y = theme_blank(), axis.text.x = theme_blank()
    with open(r_file, 'w') as f:
        f.write('library(ggplot2)\n')
        f.write("data<-read.delim('" + input_file 
                + "',header=TRUE,sep='\\t')\n")
        f.write("png('" + plot_file + "')\n")
        f.write("ggplot(data) + aes(x=Pos,y=LogRatio) + geom_point() + facet_grid(Chr~.) + opts(legend.position='none', title='Log Ratio')\n")
        f.write('dev.off()\n')
        
    os.system('R --vanilla < ' + r_file)

def plot_ref_allele_freq(snps, plot_file):
    """Use R ggplot2 to plot the reference allele frequencies for each chr"""

    input_file = 'input' + str(random.randint(0,1000))
    rtmp = 'rtmp' + str(random.randint(0,1000))

    mk_r_input(snps, input_file)
    mk_r_file(input_file, rtmp, plot_file)

def plot_allele_ratio(snps, plot_file):
    """Use R ggplot2 to plot the reference allele / min allele frequency for each chr"""

    input_file = 'input' + str(random.randint(0,1000))
    rtmp = 'rtmp' + str(random.randint(0,1000))

    mk_r_ratio_input(snps, input_file)
    mk_r_ratio_file(input_file, rtmp, plot_file.replace('png', 'ratio.png'))

    os.system('rm ' + input_file + ' ' + rtmp)                               

def get_ref_freq(str_length, call_str):
    """Return referene allele frequency at this position.
       . and , indicate a reference call"""

    ref_count = 0
    for s in call_str:
        if s in ['.', ',']:
            ref_count += 1
    return float(ref_count)/float(str_length)

def add_pos(coverage, quality, chr, pos, str_call, snps, sample):
    """Add the alleles at this chr pos to snps if the coverage and quality is high enough"""

    if coverage >= 8 and '_' not in chr and 'M' not in chr and quality > float(100):
        ref_freq = get_ref_freq(coverage, str_call)
        snps[sample][chr+':'+pos] = [ref_freq, float(1)-ref_freq]

def get_ref_nonRef_allele_freqs(afile):
    """Parse the raw exome data file to get the reference and nonRef allele frequencies.
       There are duplicate positions in the file, so I store the parsed data in a {}"""

    snps = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            for sample in samples:
                chr = sp[2]
                pos = sp[3]
                quality = float(sp[idx+1])
                coverage = int(sp[idx+2])
                call_str = sp[idx+3]
                add_pos(coverage, quality, chr, pos,
                        call_str, snps, sample)
                idx += 5
    return snps

def get_coverage(afile):
    """Print the coverage at each SNP"""

    coverages = defaultdict(dict)
    with open(afile) as f:
        for line in f:
            sp = line.split('\t')
            samples = sp[8].split('-')
            idx = 17
            chr = sp[2]
            pos = sp[3]
            for sample in samples:
                if '_' not in chr and 'M' not in chr:
                    coverage = int(sp[idx+2])
                    coverages[sample][chr+':'+pos] = int(coverage)
                idx += 5
    return coverages

data_dir = 'data/exome'
plot_dir = 'plots/cnv/'
cnv_seq_dir = 'working/cnv_seq/'
os.system('mkdir -p ' + plot_dir)
paired_data = []
for cancer, normal in global_settings.pairs:
    paired_data.append(cancer)
    paired_data.append(normal)

for afile in os.listdir(data_dir):
    if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
        file = os.path.join(data_dir, afile)
        allele_freqs = get_ref_nonRef_allele_freqs(file)
        coverages = get_coverage(file)
        plots = {}
#        for sample in paired_data:
#            plot_file = os.path.join(plot_dir, afile + '.' + sample + '.png')
#            plot_ref_allele_freq(allele_freqs[sample], plot_file)
#            plot_allele_ratio(allele_freqs[sample], plot_file)
#            file_name = os.path.join(plot_dir, afile + '.' + sample + '.sum.png')
#            os.system('montage -geometry 500 -quality 100 '
#                      + plot_file + ' '
#                      + plot_file.replace('png', 'ratio.png') + ' '
#                      + file_name)
#            plots[sample] = file_name
        for cancer, normal in global_settings.pairs:
            cov_ratio_file = os.path.join(plot_dir, afile + '.' + cancer + '.coverage.png')
            plot_coverage_ratio(coverages, cancer, normal, cov_ratio_file)
            cnv_seq_cancer_file =  os.path.join(cnv_seq_dir, afile + '.' + cancer + '.coverage.hits')
            cnv_seq_normal_file =  os.path.join(cnv_seq_dir, afile + '.' + normal + '.coverage.hits')
            mk_cvn_seq(coverages, cancer, normal, cnv_seq_cancer_file, cnv_seq_normal_file)
 #           os.system('montage -geometry 500 -quality 100 '
 #                     + plots[normal] + ' '
 #                     + plots[cancer] + ' '
 #                     + os.path.join(plot_dir, afile + '.' + cancer + '.cmp.png'))

