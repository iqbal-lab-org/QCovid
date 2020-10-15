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
for line in open('nCoV-2019.tsv'):
    name, pool, seq, *rest = line.strip().split('\t')
    if name == 'name':
        # header
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
        amplicons[name][1] = ref.index(seq) + len(seq)
    elif direction == 'LEFT':
        amplicons[name][0] = ref.index(seq)

for a in amplicons:
    start, end = amplicons[a]
    print(f"{a}\t{start}\t{end}")
