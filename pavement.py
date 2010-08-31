from paver.easy import *
import os

@task
def plt():
    """Parse data and mk LOH plots"""

    sh('python plt.py')
    sh('./liftOver working/BED_pos_hg19 data/ucsc/hg19ToHg18.over.chain working/BED_pos_hg18 working/liftOver_unmapped')
    # sh('python r_plot.py')
    
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
#             sh('montage -geometry 1000 -tile 1x1 -quality 100 '
#                + os.path.join('plots', dir, normal + '.png') + ' '
#                + os.path.join('plots', dir, cancer + '.png') + ' '
#                + os.path.join('plots', dir, normal + '_both.png'))
#             sh('rm ' + os.path.join('plots', dir, normal + '.png'))
#             sh('rm ' + os.path.join('plots', dir, cancer + '.png'))
