"""MU2A class"""

class MU2A:
    """Simple class to hold MU2A parsed data"""

    def __init__(self, mu2a_file):
        """Parses MU2A output"""

        self.data = {}
        with open(mu2a_file) as f:
            fields = f.readline().strip().split('\t')
            for line in f:
                sp = line.split('\t')
                chr = sp[2]
                pos = sp[3]
                data_dict = {}
                for idx in xrange(2, len(fields)):
                    data_dict[fields[idx]] = sp[idx].strip()
                    if data_dict[fields[idx]] == '':
                        del data_dict[fields[idx]]
                if chr not in self.data:
                    self.data[chr] = {}
                if pos not in self.data[chr]:
                    self.data[chr][pos] = []
                self.data[chr][pos].append(data_dict)

    def get_kinase_activity_mutations(self):
        """Which mutations are in a kinase?"""

        chrpos = {}
        for chr in self.data:
            for pos in self.data[chr]:
                for data_dict in self.data[chr][pos]:
                    if 'go_terms' in data_dict:
                        if 'kinase activity' in data_dict['go_terms']:
                            chrpos[chr + ':' + pos] = True
        return chrpos
                        
    def get_snp_positions(self):
        """Which mutations are already known SNPs?"""

        chrpos = {}
        for chr in self.data:
            for pos in self.data[chr]:
                for data_dict in self.data[chr][pos]:
                    if 'snpName' in data_dict:
                        chrpos[chr + ':' + pos] = True
        return chrpos
