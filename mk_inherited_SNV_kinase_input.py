chrpos = {}
with open('working/mu2a/inherited.mu2a_output.kinase_input') as f:
    for line in f:
        chr, pos, m1, m2 = line.strip().split('\t')
        chrpos[chr+':'+pos] = True
with open('working/novel_mutations/inherited.novel') as f:
    for line in f:
        chr, pos, m1, m2 = line.strip().split('\t')
        if chr+':'+pos in chrpos:
            print line.strip()
