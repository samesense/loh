"""Now I have converted hg19 to hg18, and I need to
   take the data from hg19 and assign it to points in hg18"""
import os, sys

def transform_chr(old_chr):
    """Take chrs in UCSC format, and put in quantsmooth format"""

    chr = old_chr.split('chr')[-1]
    if 'XY' == chr:
        chr = '100'
    elif 'X' == chr:
        chr = '98'
    elif 'Y' == chr:
        chr = '99'
    return chr        

def get_unmapped():
    """Get coords that I couldnt transform to hg18"""

    unmapped = {}
    with open('working/liftOver_unmapped') as f:
        for line in f:
            if line[0] != '#':
                chr, pos, junk = line.split('\t')
                unmapped[chr + ':' + pos] = True
    return unmapped

def mk_map():
    """Map hg19 to hg18. Redo chr labels"""

    unmapped = get_unmapped()
    
    map = {}
    with open('working/BED_pos_hg19') as hg19:
        with open('working/BED_pos_hg18') as hg18:
            line19 = hg19.readline()
            line18 = hg18.readline()
            while line19 != '':
                if line18.strip():
                    chr19, pos19, junk = line19.split('\t')
                    if chr19 + ':' + pos19 not in unmapped:
                        chr18, pos18, junk = line18.split('\t')
                        chr19 = transform_chr(chr19)
                        chr18 = transform_chr(chr18)
                        if chr19 != chr18:
                            sys.stderr.write('conversion problem '
                                             + chr19 + ' ' + chr18 + '\n')
                        map[chr19 + ':' + pos19] = pos18
                line19 = hg19.readline()
                line18 = hg18.readline()
    return map

def rewrite_file(file, new_file, map):
    """Redo value file for hg19 with hg18 coords"""

    with open(file) as of:
        with open(new_file, 'w') as nf:
            nf.write(of.readline())
            for line in of:
                chr, pos, val = line.split('\t')
                key = chr + ':' + pos
                # some locations were not mapped
                if key in map:
                    new_pos = map[key]
                    nf.write('%s\t%s\t%s' %
                             (chr, new_pos, val))
map = mk_map()
working_dir = 'working/'
for dir in os.listdir(working_dir):
    if 'exome' in dir:
        exome_path = os.path.join(working_dir,
                                  dir)
        for file in os.listdir(exome_path):
            if 'hg19' in file:
                oldfile = os.path.join(exome_path, file)
                newfile = os.path.join(exome_path, 
                                       file.replace('19', '18'))
                rewrite_file(oldfile, 
                             newfile,
                             map)


