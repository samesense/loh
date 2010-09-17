from paver.easy import *
import os, sys
sys.path.append('./')
import global_settings, utils

@task
def pure_aa_chg():
    """Restrict mutations to only those in aa_chg"""

    prev = 'working/exome.aa_chg.vars'
    counter = 0
    for exome in global_settings.exome_types:
        if exome != 'exome.aa_chg':
            outfile = 'working/aa_chg.pure' + str(counter) + '.vars'
            sh('comm -23 <(sort ' + prev + ') <(sort working/' + exome + '.vars) > '
               + outfile)
            prev = outfile
            counter += 1

@task
def murim_pure_aa_chg():
    """Restrict mutations to only those in aa_chg"""

    prev = 'working/murim.exome.aa_chg.vars'
    counter = 0
    for exome in global_settings.exome_types:
        if exome != 'exome.aa_chg':
            outfile = 'working/murim.aa_chg.pure' + str(counter) + '.vars'
            sh('comm -23 <(sort ' + prev + ') <(sort working/murim.' + exome + '.vars) > '
               + outfile)
            prev = outfile
            counter += 1

@task 
def mixture():
    """Find mixture between cancer and normal. Call cnv_seq before this."""

    snp_evidence = '3'
    for subdir in ('exome', 'all_non_ref_hg19'):
        sh('mkdir -p working/mixture/')
        sh('mkdir -p working/mixture/' + subdir)
        for exome_type in global_settings.exome_types:
            for cancer, normal in global_settings.pairs:
                cnv_file = 'working/cnv_seq/CNV/' + subdir + '/' + exome_type + '.' + cancer + '.cnvs'
                if os.path.exists(cnv_file):
                    sh('python find_mixture.py '
                       + 'working/murim_plot/' + subdir + '/' + exome_type.replace('.', '_') + '/hg19_murim.' + cancer + ' '
                       + cnv_file 
                       + ' > working/mixture/' + subdir + '/' + exome_type + '.' + cancer + '.mix')
        sh('python plot_mixture.py '
           + subdir + ' ' + snp_evidence)

@task
def cnv_seq():
    """Run cnv_seq"""

    sh('python cnv_plot.py')

    for subdir in ('exome', 'all_non_ref_hg19'):
        sh('mkdir -p ' + os.path.join('plots', 'cnv-seq', subdir))
        for exome in global_settings.exome_types:
            for cancer, normal in global_settings.pairs:
                cancer_hits = 'working/cnv_seq/' + subdir + '/' + exome + '.'  + cancer + '.coverage.hits'
                normal_hits = 'working/cnv_seq/' + subdir + '/' + exome + '.' + normal + '.coverage.hits'
                if utils.check_input(cancer_hits) and utils.check_input(normal_hits):
                    sh('perl cnv-seq.pl --test ' + cancer_hits + ' --ref '
                       + normal_hits + ' --genome human --log2 0.6 --p 0.001 --bigger-window 1.5 --annotate -minimum-windows 4')# --global-normalization')
                    # call loh_full to get murim data
                    sh('python cnv_seq_plot.py '
                       + exome + '.' + cancer + '.coverage.hits-vs-' + exome + '.'
                       + normal + '.coverage.hits.log2-0.6.pvalue-0.001.minw-4.cnv '
                       + 'working/murim_plot/' + subdir + '/' + exome.replace('.','_') + '/hg19_murim.' + cancer + ' '
                       + subdir + ' '
                       + 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer + '.png')
                    sh('montage -geometry 500 -quality 100 -tile 1x2 '
                       + 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer + '.png '
                       + 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer + '.murim.png '
                       + 'plots/cnv-seq/' + subdir + '/' + exome + '.' + cancer + '.both.png')

@task
def loh_homo():
    """Make Murim LOH plot for only het>homo"""

    sh('python loh.py 100 > working/loh_percents_100')
    sh('python murim_plot.py homo working/loh_mutations > working/loh_problems')
    #sh('python summary_plot.py hg19_murim plots/murim_homo/')

@task
def loh_full():
    """Make Murim LOH plot for het>any. Call this for mixture work."""

    # this makes hg19_murim stuff for mixture
    sh('python murim_plot.py all working/loh_mutations > working/loh_problems')
    sh('python summary_plot.py hg19_murim plots/murim_all/')

@task
def het():
    """Make Murim LOH plot for het->het mutations"""

    sh('python same_het.py 100 > working/het2het_percents_100')
    sh('python murim_plot.py homo working/het2het_calls > working/loh_het_problems')
    sh('python summary_plot.py hg19_murim plots/murim_het2het/')

@task
def plt():
    """Parse data and mk LOH plots"""

    sh('python get_het_blood_snps.py > working/het_snps')
    sh('python plt.py')
    sh('python summary_plot.py')
    #sh('./liftOver working/BED_pos_hg19 data/ucsc/hg19ToHg18.over.chain working/BED_pos_hg18 working/liftOver_unmapped_hg19tohg18')
    #sh('./liftOver working/BED_pos_hg18 data/ucsc/hg18ToHg17.over.chain working/BED_pos_hg17 working/liftOver_unmapped_hg18tohg17')
    #sh('python convert_BED.py hg19 hg18')
    #sh('python convert_BED.py hg18 hg17')
#    sh('python r_plot.py')
    
#     # pair up normal and cancer
#     pairs = (('yuakerN', 'yuaker'),
#              ('yuiriN', 'yuiri'),
#              ('yusanN', 'yusan'),
#              ('yusikPLX', 'yusik'),
#              ('yunoca091283P', 'yunoca091225T'))
#     dirs = ('exome_aa_chg', 'exome_intron', 
#             'exome_no_aa_chg', 'exome_unknown',
#             'exome_UTR')
#     for dir in dirs:
#         for normal, cancer in pairs:
#             sh('montage -geometry 1000 -quality 100 '
#                + os.path.join('plots', dir, normal + '.png') + ' '
#                + os.path.join('plots', dir, cancer + '.png') + ' '
#                + os.path.join('plots', dir, normal + '_both.png'))
#             sh('rm ' + os.path.join('plots', dir, normal + '.png'))
#             sh('rm ' + os.path.join('plots', dir, cancer + '.png'))
