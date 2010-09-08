"""R commands for cvn-seq plot"""
import sys, os, random

random.seed()

rtmp = 'rtmp' + str(random.randint(0,1000))
input = sys.argv[1]
murim_input = sys.argv[2]
png = sys.argv[3]


# CNV plot
with open(rtmp, 'w') as f:
    f.write("source('funcs.R')\n")
    f.write("data<-read.delim('"
            + input + "')\n")
    f.write("png('" + png + "')\n")
    f.write("plot.cnv.all(data,colour=9)\n")
    f.write('dev.off()\n')
    f.write('q()\n')

os.system('R --vanilla < ' + rtmp)


# Murim's plot
with open(rtmp, 'w') as f:
    f.write("source('funcs.R')\n")
    f.write("data<-read.delim('"
            + murim_input + "')\n")
    f.write("png('" + png.replace('png', 'murim.png') + "')\n")
    f.write("plot.murim(data,colour=9)\n")
    f.write('dev.off()\n')
    f.write('q()\n')

os.system('R --vanilla < ' + rtmp)

os.system('rm ' + rtmp)
os.system('mv ' + input + ' working/cnv_seq/CNV/')
os.system('mv ' + input.replace('cnv', 'count') + ' working/cnv_seq/CNV/')
