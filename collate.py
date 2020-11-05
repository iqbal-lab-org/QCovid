import glob
import sys

amplicons = []
samples = {}

for fn in glob.glob(f'{sys.argv[1]}/*.{sys.argv[2]}.tsv'):
    sn = fn.split('.')[0].split('/')[-1]

    d = {}
    for ln in open(fn):
        a, c = ln.strip().split()
        d[a] = c
        if a not in amplicons:
            amplicons.append(a)
    samples[sn] = d

print(','.join(['#sample',] + amplicons))

for sn in samples:
    row = list([samples[sn][a] for a in amplicons])
    print(','.join([sn,] + row))
