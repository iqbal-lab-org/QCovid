import mappy as mp
import sys

from collections import namedtuple, defaultdict
import argparse

Primer = namedtuple(
    "Primer", ["amplicon", "name", "left", "forward", "rc", "pos", "length"]
)
Read = namedtuple("Read", ["name", "seq", "qual", "comment"])
Matched = namedtuple("Matched", ["r1", "p1", "r2", "p2"])


class Stats:
    def __init__(self):
        self.total = 0
        self.neither_match = 0
        self.r2_only = defaultdict(int)
        self.r1_only = defaultdict(int)

        self.bad_ends = defaultdict(int)

        self.amplicons = defaultdict(int)
        self.matches = 0
        self.mismatch = defaultdict(int)
        self.mismatch_count = 0


class Primers:
    def __init__(self, fn):
        self.name = fn
        self.min = 100
        self.max = 0
        self.seqs = {}
        _seqs = {}

        for l in open(fn):
            amplicon, name, seq, left, forward, pos = l.strip().split("\t")
            pos = int(pos)
            forward = forward.lower() in ["t", "+", "forward", "true"]
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


def match_multi(fq1, fq2, primersets):
    for r1, r2 in zip(
        mp.fastx_read(fq1, read_comment=True), mp.fastx_read(fq2, read_comment=True)
    ):
        r1 = Read(*r1)
        r2 = Read(*r2)
        matches = {}
        for pset in primersets:
            matches[pset.name] = Matched(r1, pset.match(r1.seq), r2, pset.match(r2.seq))
        yield r1, r2, matches


def main(vargs):
    primersets = []
    for pset in vargs.primers.split(","):
        primersets.append(Primers(pset))
    aligner = mp.Aligner(vargs.ref)
    stats = defaultdict(Stats)

    for r1, r2, matches in match_multi(vargs.R1, vargs.R2, primersets):
        for pset in matches:
            m = matches[pset]
            stats[pset].total += 1
            if not m.p1 and not m.p2:
                stats[pset].neither_match += 1
                continue

            if not m.p1:
                stats[pset].r2_only[m.p2.name] += 1
                continue

            if not m.p2:
                stats[pset].r1_only[m.p1.name] += 1
                continue

            combo = "-".join(list(sorted([m.p1.name, m.p2.name])))

            if m.p1.left == m.p2.left:
                # both primers are from the same side of the template!
                stats[pset].bad_ends[combo] += 1
                continue

            if m.p1.amplicon == m.p2.amplicon:
                # good amplicon match
                stats[pset].amplicons[m.p1.amplicon] += 1
                stats[pset].matches += 1
                continue
            else:
                # mispriming
                stats[pset].mismatch[combo] += 1
                stats[pset].mismatch_count += 1

    for pset in stats:
        print(pset, stats[pset].matches, stats[pset].mismatch_count, stats[pset].total)


def match_remainder(match, aligner):
    if False:
        if match.p1 == None:
            # print("\t", match.r1.name, match.r1.seq)
            for hit in aligner.map(match.r1.seq):
                print(
                    "\t{}\t{}\t{}\t{}".format(
                        hit.ctg, hit.r_st, hit.r_en, hit.cigar_str
                    )
                )
        else:
            # print(match.r1.name, match.r1.seq)
            pass


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
