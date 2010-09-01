"""Now I have converted hg19 to hg18, and I need to
   take the data from hg19 and assign it to points in hg18.
   Bases on one chromosome can be mapped to bases on another."""
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

def get_unmapped(hg_from, hg_to):
    """Get coords that I couldnt transform to hg18"""

    unmapped = {}
    with open('working/liftOver_unmapped_' + hg_from + 'to' + hg_to) as f:
        for line in f:
            if line[0] != '#':
                chr, pos, junk = line.split('\t')
                unmapped[chr + ':' + pos] = True
    return unmapped

def mk_map(hg_from, hg_to):
    """Map hg19 to hg18. Redo chr labels"""

    unmapped = get_unmapped(hg_from, hg_to)
    
    map = {}
    with open('working/BED_pos_' + hg_from) as hg19:
        with open('working/BED_pos_' + hg_to) as hg18:
            line19 = hg19.readline()
            line18 = hg18.readline()
            while line19 != '':
                chr19, pos19, junk = line19.split('\t')
                if chr19 + ':' + pos19 not in unmapped:
                    chr18, pos18, junk = line18.split('\t')
                    if '_' not in chr18:
                        chr19 = transform_chr(chr19)
                        chr18 = transform_chr(chr18)
                        # if chr19 != chr18:
#                             sys.stderr.write('conversion problem '
#                                              + line19 + line18)
#                             sys.exit(-1) 
                        map[chr19 + ':' + pos19] = chr18 + ':' + pos18
                    line18 = hg18.readline()
                line19 = hg19.readline()
                
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
                    new_chr, new_pos = map[key].split(':')
                    nf.write('%s\t%s\t%s' %
                             (new_chr, new_pos, val))

def main(hg_from, hg_to):
    """Call to convert hgXX to hgYY"""

    map = mk_map(hg_from, hg_to)
    working_dir = 'working/'
    for dir in os.listdir(working_dir):
        if 'exome' in dir:
            exome_path = os.path.join(working_dir,
                                      dir)
            for file in os.listdir(exome_path):
                if hg_from in file:
                    oldfile = os.path.join(exome_path, file)
                    newfile = os.path.join(exome_path, 
                                           file.replace(hg_from, hg_to))
                    rewrite_file(oldfile, 
                                 newfile,
                                 map)

# enter hg_from hg_to
main(sys.argv[1], sys.argv[2])


