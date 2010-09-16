"""Use R ggplot2 to plot the mixture for all exome types and cancer/normal pairs"""
import global_settings, random, os, sys

random.seed()

subdir = sys.argv[1] # 'all_non_ref_hg19' | 'exome'
snp_cut = int(sys.argv[2])
rinput = 'rinput' + str(random.randint(0,1000))
working_dir = os.path.join('working/mixture/', subdir + '/')
with open(rinput, 'w') as f:
    f.write('Sample\tExome\tChr\tMixture\n')
    for exome_type in global_settings.exome_types:
        for cancer, normal in global_settings.pairs:
            mix_file_name = working_dir + exome_type + '.' + cancer + '.mix'
            if os.path.exists(mix_file_name):
                with open(mix_file_name) as mixfile:
                    for line in mixfile:
                        chr, snp_count, avg_diff = line.strip().split('\t')
                        if int(snp_count) > snp_cut:
                            f.write('%s\t%s\t%s\t%f\n' %
                                    (cancer, exome_type.split('.')[1], chr, 
                                     float(0.5)-float(avg_diff)))

# make R calls
rtmp = 'rtmp' + str(random.randint(0,1000))
with open(rtmp, 'w') as f:
    f.write('library(ggplot2)\n')
    f.write("data<-read.delim('" + rinput + "',header=TRUE,sep='\\t')\n")
    f.write("png('plots/" + subdir + ".mixture.png')\n")
    f.write('ggplot(data) + aes(x=Chr,y=Mixture) + geom_point() + facet_grid(Sample~Exome)\n')
    f.write('dev.off()\n')
    f.write('q()\n')

os.system('R CMD BATCH --vanilla ' + rtmp + ' tmpLog')
os.system('rm ' + rinput + ' ' + rtmp)
        

