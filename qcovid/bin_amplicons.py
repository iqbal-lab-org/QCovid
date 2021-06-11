#!/usr/bin/env python3
"""bin reads by amplicon coordinates
"""

import sys
import pysam
import argparse

# loosen the coordinate interval to catch non-softclipped adapter sequences
PADDING = 5


def load_amplicons(amplicon_bed):
    amplicons = {}
    for line in open(amplicon_bed):
        if line[0] == "#":
            continue
        fname, start, end = line.strip().split()
        amplicons[fname] = [int(start), int(end)]
    return amplicons


def bin_amplicons_se(amplicons, bam, filtered=None, reference=None):
    histogram = {}
    total = 0
    f1r2 = 0
    f2r1 = 0
    other = 0
    print(
        "\t".join(
            [
                "#amplicon",
                "f1r2",
                "f2r1",
                "r1f2",
                "r2f1",
                "total",
                "amp50",
                "amp75",
                "amp90",
            ]
        )
    )

    for amplicon in amplicons:
        histogram[amplicon] = (0, 0, 0, 0, 0, 0, 0, 0)
    for amplicon in amplicons:
        start, end = amplicons[amplicon]

        f1r2, f2r1, r1f2, r2f1, amp50, amp75, amp90, total = histogram[amplicon]
        amplicon_length = end - start

        reads = bam.fetch(reference, start, end)

        for read in reads:
            amp_cov = abs(read.template_length) / amplicon_length

            read_end = read.template_length + read.reference_start
            if read.reference_start >= start and read_end <= end:
                if amp_cov >= 0.5:
                    amp50 += 1
                if amp_cov >= 0.75:
                    amp75 += 1
                if amp_cov >= 0.9:
                    amp90 += 1
                total += 1

                if not read.is_reverse:
                    f1r2 += 1
                elif read.is_reverse:
                    r1f2 += 1  # r1' ===> <=== r2

        histogram[amplicon] = (f1r2, f2r1, r1f2, r2f1, total, amp50, amp75, amp90)
        print(
            "\t".join(
                list(
                    map(
                        str,
                        [amplicon, f1r2, f2r1, r1f2, r2f1, total, amp50, amp75, amp90],
                    )
                )
            )
        )
    return histogram


def bin_amplicons(amplicons, bam, filtered=None, reference=None):
    histogram = {}
    total = 0
    f1r2 = 0
    f2r1 = 0
    other = 0
    print(
        "\t".join(
            [
                "#amplicon",
                "f1r2",
                "f2r1",
                "r1f2",
                "r2f1",
                "total",
                "amp50",
                "amp75",
                "amp90",
            ]
        )
    )

    for amplicon in amplicons:
        histogram[amplicon] = (0, 0, 0, 0, 0, 0, 0, 0)
    for amplicon in amplicons:
        start, end = amplicons[amplicon]

        f1r2, f2r1, r1f2, r2f1, amp50, amp75, amp90, total = histogram[amplicon]
        amplicon_length = end - start

        reads = bam.fetch(reference, start, end)

        for read in reads:

            if not read.is_proper_pair:  # both mates are mapped
                continue
            if not read.is_read1:
                continue

            amp_cov = abs(read.template_length) / amplicon_length

            # branch on template orientation, r1 starts before r2
            if read.reference_start < read.next_reference_start:
                mate_end = read.template_length + read.reference_start
                if read.reference_start >= start and mate_end <= end:
                    if amp_cov >= 0.5:
                        amp50 += 1
                    if amp_cov >= 0.75:
                        amp75 += 1
                    if amp_cov >= 0.9:
                        amp90 += 1
                    total += 1

                    if (not read.is_reverse) and read.mate_is_reverse:
                        f1r2 += 1  # r1 ===> <=== r2'
                    elif read.is_reverse and (not read.mate_is_reverse):
                        r1f2 += 1  # r1' ===> <=== r2
            else:
                if read.reference_end <= end and read.next_reference_start >= start:
                    if amp_cov >= 0.5:
                        amp50 += 1
                    if amp_cov >= 0.75:
                        amp75 += 1
                    if amp_cov >= 0.9:
                        amp90 += 1
                    total += 1

                    if (not read.is_reverse) and read.mate_is_reverse:
                        r2f1 += 1  # r2' ===> <=== r1
                    elif read.is_reverse and (not read.mate_is_reverse):
                        f2r1 += 1  # r2 ===> <=== r1'

        histogram[amplicon] = (f1r2, f2r1, r1f2, r2f1, total, amp50, amp75, amp90)
        print(
            "\t".join(
                list(
                    map(
                        str,
                        [amplicon, f1r2, f2r1, r1f2, r2f1, total, amp50, amp75, amp90],
                    )
                )
            )
        )
    return histogram


def write_counts(amps_fd, histogram):
    for a in histogram:
        print("\t".join([a, str(histogram[a])]), file=amps_fd)
    amps_fd.close()


parser = argparse.ArgumentParser(
    description="Extract amplicon performance stats from name sorted bam"
)
parser.add_argument("amplicon_bed", help="bed file of amplicon positions")
# parser.add_argument('reference', help='fasta reference in same coordinate system as amplicon bed')
parser.add_argument("bam", help="name sorted bam file input")
parser.add_argument(
    "--se",
    help="input reads are single-end (Nanopore)",
    action="store_true",
    default=False,
)
parser.add_argument("--filter", help="write out new bam file with annotated reads")
parser.add_argument("--mask", help="write out amplicons which fail QC threshold")
parser.add_argument(
    "--min_coverage",
    help="minimum number of read pairs covering an amplicon to pass qc",
    default=25,
    type=int,
)
parser.add_argument(
    "--min_template_match_75",
    help="percentage of fragments which much cover at least 75% of the template",
    default=0.0,
    type=float,
)
parser.add_argument(
    "-r", "--reference", help="name of contig that reads were mapped to", default=None
)

args = parser.parse_args()


def main():
    amplicon_set = args.amplicon_bed.split(".bed")[0].split("/")[-1]

    bam = pysam.AlignmentFile(args.bam, "rb")

    filtered = None
    if args.filter:
        filtered = pysam.AlignmentFile(args.filter, "wb", template=bam)

    amplicons = load_amplicons(args.amplicon_bed)

    if args.se:
        histogram = bin_amplicons_se(
            amplicons, bam, filtered=filtered, reference=args.reference
        )
    else:
        histogram = bin_amplicons(
            amplicons, bam, filtered=filtered, reference=args.reference
        )

    if args.mask:
        bad_amps = open(args.mask, "w")
        c = 0
        t = 0.0
        for amplicon in histogram:
            _, _, _, _, x, _, _, _ = histogram[amplicon]
            t += x
            c += 1

        avg = 0.0
        if c != 0:
            avg = t / c

        for amplicon in histogram:
            _, _, _, _, coverage, _, a75, _ = histogram[amplicon]
            # if x == 0 or x < (avg / 4) or (a50 / x) < 0.5:
            if coverage < args.min_coverage or a75 < args.min_template_match_75 - 75:
                print(amplicon, file=bad_amps)
        bad_amps.close()


if __name__ == "__main__":
    main()
