import mappy as mp
import sys


class Primers:
    def __init__(self, fn):
        self.min = 100
        self.max = 0
        self.seqs = {}
        _seqs = {}

        for l in open(fn):
            amplicon, name, seq, left, forward, pos = l.strip().split(",")
            pos = int(pos)
            forward = forward.lower() in ["t", "true"]
            left = left.lower() in ["left", "true", "t"]

            length = len(seq)
            if length > self.max:
                self.max = length
            if length < self.min:
                self.min = length

            _seqs[seq] = (amplicon, name, left, forward, True, pos, length)
            _seqs[mp.revcomp(seq)] = (amplicon, name, left, forward, False, pos, length)

        for k, v in _seqs.items():
            self.seqs[k[: self.min]] = v

    def match(self, seq):
        if seq[: self.min] in self.seqs:
            return self.seqs[seq[: self.min]]
        else:
            return None


def readpairs(fq1, fq2, matchfn):
    for r1, r2 in zip(mp.fastx_read(fq1), mp.fastx_read(fq2)):
        yield (*r1, matchfn(r1[1]), (*r2, matchfn(r2[1])))


primers = Primers(sys.argv[1])

for r1 in readpairs(sys.argv[2], sys.argv[3], primers.match):
    print(r1)
