import sys
import mappy

#a = mappy.Aligner(sys.argv[1])

for line in open('midnight-1200.tsv'):
    name, seq, _, _, _, _, pos, _ = line.strip().split('\t')
    if name == "Primer Name":
        continue

    r = name.split('_')
    if len(r) == 4:
        d = r[-1]
    else:
        assert False

    amplicon = '_'.join(r[:3])

    left = d == "LEFT"
    forward = left 

    print('\t'.join([amplicon, name, seq, str(left), str(forward), str(pos)]))
