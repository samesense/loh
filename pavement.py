from paver.easy import *

@task
def plt():
    """Parse data and mk LOH plots"""

    sh('python plt.py')
    sh('python r_plot.py')
       
