"""Use R ggplot2 to plot the mixture for all exome types and cancer/normal pairs"""
import global_settings, random, os

random.seed()

rinput = 'rinput' + str(random.randint(0,1000))
working_dir = 'working/mixture/'
with open(rinput, 'w') as f:
    f.write('Sample\tExome\tChr\tMixture\n')
    for exome_type in global_settings.exome_types:
        for cancer, normal in global_settings.pairs:
            with open(working_dir 
                      + exome_type + '.' 
                      + cancer + '.mix') as mixfile:
                for line in mixfile:
                    chr, snp_count, avg_diff = line.strip().split('\t')
                    f.write('%s\t%s\t%s\t%f\n' %
                            (cancer, exome_type, chr, 
                             float(0.5)-float(avg_diff)))

# make R calls
rtmp = 'rtmp' + str(random.randint(0,1000))
with open(rtmp, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("data<-read.delim('" + rinput + "',header=TRUE,sep='\\t')\n")
    f.write("png('plots/mixture.png')\n")
    f.write('ggplot(data) + aes(x=Chr,y=Mixture) + geom_point() + facet_grid(Sample~Exome)\n')
    f.write('dev.off()\n')
    f.write('q()\n')

os.system('R --vanilla < ' + rtmp)
os.system('rm ' + rinput + ' ' + rtmp)
        
