import mappy as mp
import sys

from collections import namedtuple
import argparse

Primer = namedtuple(
    "Primer", ["amplicon", "name", "left", "forward", "rc", "pos", "length"]
)
Read = namedtuple("Read", ["name", "seq", "qual", "comment"])
Matched = namedtuple("Matched", ["r1", "p1", "r2", "p2"])


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

            _seqs[seq] = Primer(amplicon, name, left, forward, True, pos, length)
            _seqs[mp.revcomp(seq)] = Primer(
                amplicon, name, left, forward, False, pos, length
            )

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
        r1 = Read(*r1)
        r2 = Read(*r2)
        yield Matched(r1, matchfn(r1.seq), r2, matchfn(r2.seq))


def main(vargs):
    primers = Primers(vargs.primers)
    aligner = mp.Aligner(vargs.ref)

    for match in readpairs(vargs.R1, vargs.R2, primers.match):
        if match.p1 == None:
            print("\t", match.r1.name, match.r1.seq)
            # for hit in aligner.map(r1.seq):

            #    print("\t{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
        else:
            print(match.r1.name, match.r1.seq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="identify amplicons by primer sequences and alignment positions"
    )
    parser.add_argument("ref")
    parser.add_argument("primers")
    parser.add_argument("R1")
    parser.add_argument("R2")
    args = parser.parse_args()
    main(args)
