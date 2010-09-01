"""This plot groups all chromosomes together.
   This is a remake of Murim's plot."""
import sys, os, random
from collections import defaultdict

random.seed()

def mk_input(plt_output, title):
    """Make input for a full chromosome plot.
       Make new mutation positions, and put
       lines for chormosomes."""

    rfile = 'tmpr' + title + str(random.randint(0,1000))
    rlines = 'tmprlines' + title + str(random.randint(0,1000))
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

def mk_r(points_input, lines_input, rfile, plot_file, title):
    """Add code to rfile to make a plot"""
    
    rfile.write("png('" + plot_file + "')\n")
    rfile.write("points <- read.delim('"
                + points_input + "',header=FALSE,sep='\\t')\n")
    rfile.write("plot(points,pch='o',col='red',xaxt='n',xlab='',ylab='',main='" + title + "')\n")
    rfile.write("v_line <- read.delim('"
                + lines_input + "',header=FALSE,sep='\\t')\n") 
    rfile.write("for (idx in 1:dim(v_line)[1]) {abline(v=v_line[idx,],lty=2)}\n")
    rfile.write('dev.off()\n')


file_type = sys.argv[1] # hg19 or hg19_murim

working_dir = 'working'
plot_dir = 'plots'
rinput = 'rtmp'
rm_ls = [rinput]
with open(rinput, 'w') as rout:
    for dir in os.listdir(working_dir):
        if 'exome' in dir:
            subdir = os.path.join(working_dir, dir)
            new_plot_dir = os.path.join(plot_dir, dir)
            os.system('mkdir -p ' + new_plot_dir)
            for filename in os.listdir(subdir):
                if file_type in filename:
                    title = dir + ':' + filename.split('.')[-1] + ':' + file_type
                    path = os.path.join(subdir,
                                        filename)
                    rfile, rlines = mk_input(path, title.replace(':', '_'))
                    rm_ls.append(rfile)
                    rm_ls.append(rlines)
                    mk_r(rfile, rlines, rout, 
                         os.path.join(new_plot_dir,
                                      filename.split('.')[-1] 
                                      + '.' + file_type + '.summary.png'),
                         title)
os.system('R --vanilla < rtmp')
os.system('rm ' + ' '.join(rm_ls))
for dir in os.listdir(working_dir):
    if 'exome' in dir:
        new_plot_dir = os.path.join(plot_dir, dir)
        os.system('montage -geometry 500 -quality 100 '
                  + os.path.join(new_plot_dir, '*' + file_type + '.summary*') + ' '
                  + os.path.join(new_plot_dir, dir + '.' + file_type + '.ALL.png'))
