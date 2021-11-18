#!/usr/bin/env python3

import mappy as mp
import sys
import json

from collections import namedtuple, defaultdict
import argparse

from qcovid.primers import Primer, Primers, Read, Matched, print_read


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
    total = 0
    for r1, r2, matches in match_multi(vargs.R1, vargs.R2, primersets):
        total += 1
        match_mode = "unmatched"
        combo = None

        for pset in matches:
            m = matches[pset]
            stats[pset].total += 1
            if not m.p1 and not m.p2:
                stats[pset].neither_match += 1
                continue

            if not m.p1:
                stats[pset].r2_only[m.p2.name] += 1
                match_mode = "no_r1"
                continue

            if not m.p2:
                stats[pset].r1_only[m.p1.name] += 1
                match_mode = "no_r2"
                continue

            combo = "-".join(list(sorted([m.p1.name, m.p2.name])))

            if m.p1.left == m.p2.left:
                # both primers are from the same side of the template!
                stats[pset].bad_ends[combo] += 1
                match_mode = "same_side"
                continue

            if m.p1.amplicon == m.p2.amplicon:
                # good amplicon match
                stats[pset].amplicons[m.p1.amplicon] += 1
                stats[pset].matches += 1
                match_mode = "exact"
                continue
            else:
                # mispriming
                stats[pset].mismatch[combo] += 1
                stats[pset].mismatch_count += 1
                match_mode = "mismatch"

        if vargs.prefix:
            if match_mode not in ["exact"] and filter_reads:
                continue

            r1.name += f" qcovid-v0.1:{combo}:{match_mode}:{plen}"
            r2.name += f" qcovid-v0.1:{combo}:{match_mode}:{plen}"

            print_read(r1, file=r1fd)
            print_read(r2, file=r1fd)

    if total == 0:
        if args.json:
            print(json.dumps({"status": "failure", "message": "no reads in input"}))
        print("No reads in input", file=sys.stderr)
        exit(0)

    firstn = None
    secondn = None
    first = 0
    second = 0
    for pset in stats:
        print(
            pset,
            stats[pset].matches,
            stats[pset].mismatch_count,
            stats[pset].total,
            file=sys.stderr,
        )
        ms = stats[pset].matches

        if ms >= first:
            firstn = pset
            first = ms

        if ms >= second and ms < first:
            secondn = pset
            second = ms

    f = first / total
    s = second / total

    if f >= 0.7 and s <= 0.3:
        results = {"status": "success", "read_pairs": total, "primer_set": firstn}
        print(f"selected primerset {firstn}", file=sys.stderr)
        if args.json:
            print(json.dumps(results))
        exit(0)

    results = {
        "status": "failure",
        "read_pairs": total,
        "message": "ambiguous primer set",
    }
    if args.json:
        print(json.dumps(results))
    exit(0)


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
    parser.add_argument("--json", required=False, action="store_true")
    parser.add_argument("--prefix", required=False)
    parser.add_argument("ref")
    parser.add_argument("primers")
    parser.add_argument("R1")
    parser.add_argument("R2")
    args = parser.parse_args()
    main(args)
