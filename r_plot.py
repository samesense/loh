import os, sys

def do_file(filename):
    new_file = filename.split('.')[-1] + '.png'
    with open('rtmp', 'w') as f:
        f.write("png('plots/" + new_file + "')\n")
        f.write('library(quantsmooth)\n')
        f.write("MapInfo<-read.delim('working/" + filename + "',header=TRUE,sep='\t')\n")
        f.write("chrompos<-prepareGenomePlot(MapInfo[,1:2],paintCytobands=TRUE,organism='hsa')\n")
        f.write('for (i in 1:11) {segments(0,i+.5,3000000000,i+.5,col="blue",lwd=1)}\n')
        f.write("points(chrompos[,2],chrompos[,1]+MapInfo[,3],pch='.',col='red')\n")
        f.write("text(300000000/2,12.5,'" + filename.split('.')[-1] + "')\n")
        f.write('dev.off()\n')
    os.system('R --vanilla < rtmp')

for filename in os.listdir('working/'):
    if 'chr_pos' in filename:
        do_file(filename)
os.system('rm rtmp')

