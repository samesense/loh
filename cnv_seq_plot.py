"""R commands for cvn-seq plot"""
import sys, os, random

rtmp = 'rtmp' + str(random.randint(0,1000))
input = sys.argv[1]
png = sys.argv[2]

with open(rtmp, 'w') as f:
    f.write('library(cnv)\n')
    f.write("data<-read.delim('"
            + input + "')\n")
    f.write("png('" + png + "')\n")
    f.write("plot.cnv.all(data,colour=9)\n")
    f.write('dev.off()\n')
    f.write('q()\n')

os.system('R --vanilla < ' + rtmp)
os.system('rm ' + rtmp + ' ' + input + ' ' + input.replace('cnv', 'count'))

