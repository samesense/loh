"""This plot groups all chromosomes together.
   This is a remake of Murim's plot."""
import sys, os, random
from collections import defaultdict

random.seed()

def mk_input(plt_output):
    """Make input for a full chromosome plot.
       Make new mutation positions, and put
       lines for chormosomes."""

    rfile = 'tmpr' + str(random.randint(0,1000))
    rlines = 'tmprlines' + str(random.randint(0,1000))
    chr2posLs = defaultdict(list)
    with open(plt_output) as f:
        f.readline()
        for line in f:
            (chr, pos, val) = [x.strip() for x in line.split('\t')]
            chr2posLs[chr].append((int(pos), val))
    new_counter = 0
    chrs = chr2posLs.keys()
    chrs.sort()
    with open(rlines, 'w') as input_lines:
        with open(rfile, 'w') as input_outfile:
            for chr in chrs:
                input_lines.write(str(new_counter) + '\n')
                chr2posLs[chr].sort()
                for pos, val in chr2posLs[chr]:
                    input_outfile.write(str(new_counter) + '\t' + val + '\n')
                    new_counter += 1
    return (rfile, rlines)

def mk_r(points_input, lines_input, rfile, plot_file):
    """Add code to rfile to make a plot"""
    
    rfile.write("png('" + plot_file + "')\n")
    rfile.write("points <- read.delim('"
                + points_input + "',header=FALSE,sep='\\t')\n")
    rfile.write("plot(points,pch='.')\n")
    rfile.write("v_line <- read.delim('"
                + lines_input + "',header=FALSE,sep='\\t')\n") 
    rfile.write("for (idx in 1:dim(v_line)[1]) {abline(v=v_line[idx,],lty=2)}\n")
    rfile.write('dev.off()\n')

rinput = 'rtmp'
rm_ls = [rinput]
rfile, rlines = mk_input('working/exome_aa_chg/hg19_chr_pos.yusik')
rm_ls.append(rfile)
rm_ls.append(rlines)
with open(rinput, 'w') as rout:
    mk_r(rfile, rlines, rout, 'tmp_plt.png')
os.system('R --vanilla < rtmp')
os.system('rm ' + ' '.join(rm_ls))
