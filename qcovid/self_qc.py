#!/usr/bin/env python3
"""Self-consistent alignment QC.
This is a wrapper around bamtools' pileup.
"""
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(
    description="Self-consistency QC stats for covid assemblies"
)
parser.add_argument("assembly", help="target assembly")
parser.add_argument("reads", help="reads mapped to assembly")
parser.add_argument(
    "--fasta", help="output masked fasta to stdout", action="store_true"
)
parser.add_argument("--freq_threshold", default=0.5, type=float)
args = parser.parse_args()


def main():
    bam = pysam.AlignmentFile(args.reads, "rb")
    ref = pysam.FastaFile(args.assembly)

    if not args.fasta:
        print(
            "\t".join(
                [
                    "reference",
                    "position",
                    "ref",
                    "ref_depth",
                    "total_depth",
                    "n_segments",
                    "n_aligned",
                    "freq",
                    "n_segments/n_aligned",
                    "insertions",
                    "deletions",
                    "bases",
                    "r1forward",
                    "r1reverse",
                    "r2forward",
                    "r2reverse",
                    "secondary_alignments",
                    "not_paired",
                ]
            )
        )

    refseq = ""
    orig_name = ""
    for line in open(args.assembly):
        if line[0] == ">":
            orig_name = line[1:].strip()
            continue
        refseq += line.strip()
    masked = list(refseq)
    for pc in bam.pileup(stepper="samtools", fastafile=ref):
        dels = 0
        ins = 0
        starts = 0
        ends = 0
        secondaries = 0
        r1forward = 0
        r1reverse = 0
        r2reverse = 0
        r2forward = 0
        not_paired = 0

        d = {
            "*": 0,
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0,
            "+": 0,
            "-": 0,
            ",": 0,
            ".": 0,
            "N": 0,
        }
        d_total = 0

        ref_base = ref.fetch(pc.reference_name, pc.reference_pos, pc.reference_pos + 1)
        assert ref_base in d
        primer_len = 25  # hack TODO: either max of all primers or lookup
        # if ref_base not in d:
        #    ref_base = "N"
        for read in pc.pileups:
            if not read.alignment.cigartuples:
                # unmapped read
                continue

            if read.alignment.is_secondary:
                secondaries += 1
                continue

            try:
                if (not read.alignment.is_reverse) and read.alignment.is_read2:
                    primer_base = read.query_position > read.query_length - primer_len
                    r2forward += 1
                if (not read.alignment.is_reverse) and read.alignment.is_read1:
                    primer_base = read.query_position > read.query_length - primer_len
                    r1forward += 1
                if (read.alignment.is_reverse) and read.alignment.is_read2:
                    primer_base = read.query_position < primer_len
                    r2reverse += 1
                if (read.alignment.is_reverse) and read.alignment.is_read1:
                    primer_base = read.query_position < primer_len
                    r1reverse += 1
            except AttributeError as e:
                # pysam needs to be removed from this project
                # I think
                # print(f"warning: {e}", file=sys.stderr)
                continue

            if not read.alignment.is_proper_pair:
                not_paired += 1

            if read.indel < 0:
                dels += 1

            if read.indel > 0:
                ins += 1

            if read.is_del:
                pass

            elif read.is_refskip:
                pass
            else:
                base = read.alignment.query_sequence[read.query_position]
                if base not in d:
                    base = "N"

                if not primer_base:
                    d[base] += 1
                    d_total += 1

        freq = 0.0
        num_aligned = pc.get_num_aligned()
        ref_count = d[ref_base]
        if d_total == 0:
            print(
                pc.reference_pos, "is empty. entirely inside primer?", file=sys.stderr
            )
            pass
        else:
            freq = ref_count / num_aligned

            if freq >= args.freq_threshold:
                continue

        bases = []
        for k in d:
            bases.append(f"{k}:{d[k]}")
        bases = ";".join(bases)
        if not args.fasta:
            print(
                "\t".join(
                    map(
                        str,
                        [
                            pc.reference_name,
                            pc.reference_pos + 1,
                            ref_base,
                            ref_count,
                            d_total,
                            pc.nsegments,
                            num_aligned,
                            freq,
                            num_aligned / pc.nsegments,
                            ins,
                            dels,
                            bases,
                            r1forward,
                            r1reverse,
                            r2forward,
                            r2reverse,
                            secondaries,
                            not_paired,
                        ],
                    )
                )
            )
        else:
            masked[pc.reference_pos] = "N"

    if args.fasta:
        print(f">{orig_name}_masked")
        print("".join(masked))


if __name__ == "__main__":
    main()
