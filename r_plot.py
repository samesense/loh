"""Call R functions to make LOH plots"""
import os, sys, random

random.seed()

def write_file(path, filename, plot_dir, title, rfile):
    """Accumulate R commands in one file. Use for making all plots at once"""

    new_file = os.path.join(plot_dir,
                            filename.split('.')[-1] + '.png')
    rfile.write("png('" + new_file + "')\n")
    rfile.write("MapInfo<-read.delim('" + path + "',header=TRUE,sep='\\t')\n")
    rfile.write("chrompos<-prepareGenomePlot(MapInfo[,1:2],paintCytobands=TRUE,organism='hsa',sexChromosomes=TRUE)\n")
    rfile.write('for (i in 1:12) {segments(0,i+.5,3000000000,i+.5,col="blue",lwd=1)}\n')
    rfile.write("points(chrompos[,2],chrompos[,1]+MapInfo[,3],pch='.',col='red')\n")
    rfile.write("text(300000000/2,13.5,'" + title + "')\n")
    rfile.write('dev.off()\n')

def do_file(path, filename, plot_dir, title):
    """Make one plot"""

    new_file = os.path.join(plot_dir,
                            filename.split('.')[-1] + '.png')
    with open('rtmp', 'w') as f:
        f.write("png('" + new_file + "')\n")
        f.write('library(quantsmooth)\n')
        f.write("MapInfo<-read.delim('" + path + "',header=TRUE,sep='\\t')\n")
        f.write("chrompos<-prepareGenomePlot(MapInfo[,1:2],paintCytobands=TRUE,organism='hsa',sexChromosomes=TRUE)\n")
        f.write('for (i in 1:12) {segments(0,i+.5,3000000000,i+.5,col="blue",lwd=1)}\n')
        f.write("points(chrompos[,2],chrompos[,1]+MapInfo[,3],pch='.',col='red')\n")
        f.write("text(300000000/2,13.5,'" + title + "')\n")
        f.write('dev.off()\n')
    os.system('R --vanilla < rtmp')

working_dir = 'working'
plot_dir = 'plots'
rfile = 'rtmp'
with open(rfile, 'w') as rfile:
    rfile.write('library(quantsmooth)\n')
    for dir in os.listdir(working_dir):
        if 'exome' in dir:
            subdir = os.path.join(working_dir, dir)
            new_plot_dir = os.path.join(plot_dir, dir)
            os.system('mkdir -p ' + new_plot_dir)
            for filename in os.listdir(subdir):
                title = dir + ':' + filename.split('.')[-1]
                path = os.path.join(subdir,
                                    filename)
                write_file(path, filename, new_plot_dir, title, rfile)
os.system('R --vanilla < rtmp')
os.system('rm rtmp')

