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
    for r1, r2 in zip(
        mp.fastx_read(fq1, read_comment=True), mp.fastx_read(fq2, read_comment=True)
    ):
        yield (*r1, matchfn(r1[1])), (*r2, matchfn(r2[1]))


primers = Primers(sys.argv[1])
aligner = mp.Aligner(sys.argv[2])


for r1, r2 in readpairs(sys.argv[3], sys.argv[4], primers.match):
    if r1[3] == None:
        print("\t", r1[0], r1[3])
        # for hit in aligner.map(r1[1]):

        #    print("\t{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
    else:
        print(r1[0], r1[3])
