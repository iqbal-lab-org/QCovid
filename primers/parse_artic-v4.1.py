import sys
import mappy

#a = mappy.Aligner(sys.argv[1])

for line in open('V4.1_alts.primer.bed'):
    _, pos, _, name, _, strand, seq = line.strip().split('\t')

    r = name.split('_')
    if len(r) == 3:
        d = r[-1]
    elif len(r) == 4:
        d = r[-2]

    amplicon = '_'.join(r[:2])

    left = d == "LEFT"
    forward = strand == "+"
    print('\t'.join([amplicon, name, seq, str(left), str(forward), str(pos)]))
