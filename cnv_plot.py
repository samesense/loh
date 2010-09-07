"""Make the copy number plot from
   High-resolution genomic profiling
   of chromosomal abberations using
   Infinium whole-genome genotyping.
   Make a plot for normal and cancer.
"""
import os, random, sys, math
from collections import defaultdict

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

data_dir = 'data/exome'
plot_dir = 'plots/cnv/'
os.system('mkdir -p ' + plot_dir)
for afile in os.listdir(data_dir):
    if not 'sun' in afile: # ignore exome.aa_chg.sun b/c is it redundant
        file = os.path.join(data_dir, afile)
        allele_freqs = get_ref_nonRef_allele_freqs(file)
        for sample in ['yusan', 'yusanN']:
#            ref_allele_freq, min_allele_freq = allele_freqs[chr_pos]
            plot_file = os.path.join(plot_dir, afile + '.' + sample + '.png')
            plot_ref_allele_freq(allele_freqs[sample], plot_file)
            plot_allele_ratio(allele_freqs[sample], plot_file)
            os.system('montage -geometry 500 -quality 100 '
                      + plot_file + ' '
                      + plot_file.replace('png', 'ratio.png') + ' '
                      + os.path.join(plot_dir, afile + '.' + sample + '.sum.png'))
#            print('%s\t%s\t%s\t%f'
#                  % (afile, chr, pos,
#                     ref_allele_freq)) 
