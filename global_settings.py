homo_bases = ('A', 'T', 'G', 'C')

het_bases = ('S', 'Y', 'M', 'R', 'K', 'W')

# cancer, normal
pairs = (('yuaker', 'yuakerN'),
         ('yuiri', 'yuiriN'),
         ('yuiskia', 'yuiskiaN'),
         ('yunoca091225T', 'yunoca091283P'),
         ('yusan', 'yusanN'))

exome_types = ('exome.aa_chg', 'exome.intron', 
               'exome.no_aa_chg', 'exome.unknown',
               'exome.UTR')

# used to make all_non_ref data look like exome data
alias = {{'yusanT':'yusan'},
         {'yusanN':'yusanN'}}
