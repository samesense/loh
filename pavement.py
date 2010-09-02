from paver.easy import *
import os

@task
def loh():
    """Make Murim LOH plot"""

    sh('python loh.py > working/loh_percents')
    sh('python murim_plot.py')
    sh('python summary_plot.py hg19_murim')

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
