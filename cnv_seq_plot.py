"""R commands for cnv-seq plot & Murim's plot. I compare the two to find the mixture between normal and cancer cells. This also outputs the CNVs to working/cnv-seq/CNV/"""
import sys, os, random
import utils

random.seed()

rtmp = 'rtmp' + str(random.randint(0,1000))
input = sys.argv[1]
murim_input = sys.argv[2]
subdir = sys.argv[3] # exome | all_non_ref_hg19
png = sys.argv[4]

# CNV plot
os.system('mkdir -p working/cnv_seq/CNV')
os.system('mkdir -p working/cnv_seq/CNV/' 
          + subdir + '/')
os.system('rm working/cnv_seq/CNV/' 
          + subdir + '/' 
          + input.split('.coverage')[0]  + '.cnvs')
with open(rtmp, 'w') as f:
    f.write("source('funcs.R')\n")
    f.write('library(cnv)\n')
    f.write("data<-read.delim('"
            + input + "')\n")
    f.write("png('" + png + "')\n")
    f.write("plot.cnv.all.perry(data,colour=9)\n")
    f.write('dev.off()\n')
    f.write("cnv.print(data, file='working/cnv_seq/CNV/" + subdir + '/' + input.split('.coverage')[0]  + ".cnvs')\n")
    f.write('q()\n')

if utils.check_input(input):
    os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')

# Murim's plot
with open(rtmp, 'w') as f:
    f.write("source('funcs.R')\n")
    f.write("data<-read.delim('"
            + murim_input + "')\n")
    f.write("png('" + png.replace('png', 'murim.png') + "')\n")
    f.write("plot.murim(data,colour=9)\n")
    f.write('dev.off()\n')
    f.write('q()\n')

if utils.check_input(murim_input):
    os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')

os.system('rm ' + rtmp + ' tmpLog')
os.system('mv ' + input + ' working/cnv_seq/CNV/' + subdir + '/')
os.system('mv ' + input.replace('cnv', 'count') + ' working/cnv_seq/CNV/' + subdir + '/')
