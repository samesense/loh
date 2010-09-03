from paver.easy import *
import os

@task
def loh_homo():
    """Make Murim LOH plot for only het>homo"""

    sh('python loh.py 100 > working/loh_percents_100')
    sh('python murim_plot.py homo working/loh_mutations > working/loh_problems')
    sh('python summary_plot.py hg19_murim plots/murim_homo/')

@task
def loh_full():
    """Make Murim LOH plot for het>any"""

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
