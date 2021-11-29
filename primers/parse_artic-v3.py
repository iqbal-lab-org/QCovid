import sys
import mappy

#a = mappy.Aligner(sys.argv[1])

poss = {}
for line in open('artic-v3.bed'):
    _, start, end, name, _, _ = line.strip().split()
    poss[name] = (start, end)

for line in open('artic-v3.tsv'):
    name, amplicon, seq, *r = line.strip().split('\t')
    if name == "name":
        continue

    r = name.split('_')
    if len(r) == 3:
        d = r[-1]
    elif len(r) == 4:
        d = r[-2]

    amplicon = '_'.join(r[:2])
    left = d == "LEFT"
    forward = d == "LEFT"
    pos = poss[name][0]
    print('\t'.join([amplicon, name, seq, str(left), str(forward), str(pos)]))
