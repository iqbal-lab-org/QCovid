"""Preprocess ARTIC-v3 primer set for QCovid amplicon QC pipeline
this required manual intervention in the source nCoV-2019.tsv for
nCoV-2019_13_pool1 TCGCACAAATGTCTACTTAGCTGT LEFT
when using MN908947.
"""
ref = ""
ref_file = "../MN908947.fasta"
for line in open(ref_file):
    if line[0] == '>':
        continue
    ref += line.strip()
print(f"#Reference: {ref_file} ({len(ref)} bp.)")

def rc(seq):
    c = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(list(map(lambda x: c[x], reversed(seq))))

amplicons = {}
mix = None
for line in open('nCoV-2019-MN908947.tsv'):
    name, pool, seq, *rest = line.strip().split('\t')
    if name == 'name':
        # header line
        continue

    # the primer set includes alternative primers. this is not
    # relevent when we're just trying to resolve the start and end
    # positions
    if 'alt' in name:
        continue

    if pool == 'nCoV-2019_1':
        mix = 1
    elif pool == 'nCoV-2019_2':
        mix = 2
    else:
        raise ValueError

    _, number, direction = name.split('_')
    name = f"nCoV-2019_{number}_pool{mix}"

    if name not in amplicons:
        amplicons[name] = [None, None]
    
    if direction == 'RIGHT':
        seq = rc(seq)
        if seq not in ref:
            print(name, seq, direction, "not found")
            raise ValueError 
        amplicons[name][1] = ref.index(seq) + len(seq)
    elif direction == 'LEFT':
        if seq not in ref:
            print(name, seq, direction, "not found")
            raise ValueError
        amplicons[name][0] = ref.index(seq)

for a in amplicons:
    start, end = amplicons[a]
    print(f"{a}\t{start}\t{end}")
