ref = ""
ref_file = "MN908947.fasta"
for line in open(ref_file):
    if line[0] == '>':
        continue
    ref += line.strip()
print(f"#Reference: {ref_file} ({len(ref)} bp.)")

def rc(seq):
    c = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(list(map(lambda x: c[x], reversed(seq))))

mix = 1
amplicons = {}
for line in open('dutch_primal_amplicons.txt'):
    if line.strip() == 'Mix 2':
        mix = 2
        continue

    name, seq, d, *comment = line.strip().split()

    name = '_'.join(name.split('_')[:-1])
    name += f"_pool{mix}"

    if name not in amplicons:
        amplicons[name] = [None, None]

    if d == 'REV':
        seq = rc(seq)
        amplicons[name][1] = ref.index(seq) + len(seq)
    else:
        amplicons[name][0] = ref.index(seq)

#    if seq not in ref:
#        print('missing seq:', seq)
    
#    print(f"{name} {ref.index(seq)}")

for a in amplicons:
    start, end = amplicons[a]
    print(f"{a}\t{start}\t{end}")
