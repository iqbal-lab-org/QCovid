import json
import argparse

from collections import defaultdict, namedtuple
from primers import ReadData
import pysam


class Stats:
    def __init__(self):
        self.vs = defaultdict(Vstats)

    def update(self, variant, good, in_primer, reverse, second):
        self.vs[variant][(good, in_primer, reverse, second)] += 1


def main(vargs):
    bam = pysam.AlignmentFile(vargs.bam, "rb")
    ref = pysam.FastaFile(args.ref)
    positions = defaultdict(Stats)
    t = 0
    for read in pc.pileups:
        t += 1
        rd = ReadData(read.query_name)
        if rd.good:
            c += 1
        if rd.amplicon:
            c += 1
        else:
            pass

        positions[pos].update(
            variant, rd.good, in_primer, read.is_reverse, read.is_second
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="identify amplicons by primer sequences and alignment positions"
    )
    parser.add_argument("ref")
    parser.add_argument("primers")
    parser.add_argument("bam")
    args = parser.parse_args()
    main(args)
