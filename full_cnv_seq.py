"""Use pileup to run CNV-seq, if possible"""
import mixture_runner

def mk_input(infile, outfile):
    """Write CNV-seq input"""

    with open(outfile, 'w') as outf:
        with open(infile) as f:
            for line in f:
                sp = line.split('\t')
                chr = sp[0]
                pos = sp[1]
                cov = int(sp[7])
                mixture_runner.mk_cnv_seq_input_chrpos(cov, chr, pos, outf)

mk_input('data/pileup_hg19/yuakerT.pileup',
         'working/cnv_seq/yuakerT.hits')
mk_input('data/pileup_hg19/yuakerN.pileup',
         'working/cnv_seq/yuakerN.hits')
            
