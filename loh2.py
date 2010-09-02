"""Im not getting the same counts btwn loh.py and somatic_changes.py
   This parses out LOH from the somatic_changes output.
   Found the problem. I left out a C in the cancer calls"""

def is_loh(cancer, normal):
    """Return true if normal->cancer is LOH"""

    if normal in ('S', 'Y', 'M', 'R', 'K', 'W') and cancer in ('A', 'T', 'G', 'C'):
        return True
    
    return False

count = 0
with open('working/somatic_mutations') as f:
    for line in f:
        exome_type, normal, cancer, n_call, c_call = line.strip().split('\t')
        if is_loh(c_call, n_call):
            count += 1
print count
