import os, sys

def do_file(path, filename, plot_dir):
    new_file = os.path.join(plot_dir,
                            filename.split('.')[-1] + '.png')
    with open('rtmp', 'w') as f:
        f.write("png('" + new_file + "')\n")
        f.write('library(quantsmooth)\n')
        f.write("MapInfo<-read.delim('" + path + "',header=TRUE,sep='\\t')\n")
        f.write("chrompos<-prepareGenomePlot(MapInfo[,1:2],paintCytobands=TRUE,organism='hsa')\n")
        f.write('for (i in 1:11) {segments(0,i+.5,3000000000,i+.5,col="blue",lwd=1)}\n')
        f.write("points(chrompos[,2],chrompos[,1]+MapInfo[,3],pch='.',col='red')\n")
        f.write("text(300000000/2,12.5,'" + filename.split('.')[-1] + "')\n")
        f.write('dev.off()\n')
    os.system('R --vanilla < rtmp')

working_dir = 'working'
plot_dir = 'plots'
for dir in os.listdir(working_dir):
    if 'exome_aa_chg' in dir:
        subdir = os.path.join(working_dir, dir)
        new_plot_dir = os.path.join(plot_dir, dir)
        os.system('mkdir -p ' + new_plot_dir)
        for filename in os.listdir(subdir):
            path = os.path.join(subdir,
                                filename)
            do_file(path, filename, new_plot_dir)
        sys.exit(0)
os.system('rm rtmp')

