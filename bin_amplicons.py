"""bin reads by amplicon coordinates
"""

import sys
import pysam
import argparse

# loosen the coordinate interval to cast non-softclipped adapter sequences
PADDING = 5

def load_amplicons(amplicon_bed):
    amplicons = {}
    for line in open(amplicon_bed):
        if line[0] == '#':
            continue
        fname, start, end = line.strip().split()
        amplicons[fname] = [int(start), int(end)]
    return amplicons

def bin_amplicons(amplicons, bam, filtered=None):
    histogram = {}
    total = 0
    for amplicon in amplicons:
        start, end = amplicons[amplicon]

        # ensure coords are sorted
        rc = 0
        for read in bam.fetch():
            if not read.reference_start or not read.reference_end:
                continue
            if read.mate_is_unmapped:
                continue
            if read.reference_start > (start - PADDING) and read.reference_end < (end + PADDING) and read.next_reference_start > start:
                rc += 1
                if filtered:
                    filtered.write(read)
                total += 1
        histogram[amplicon] = rc
    return histogram

def write_counts(amps_fd, histogram):
    for a in histogram:
        print('\t'.join([a, str(histogram[a])]), file=amps_fd)
    amps_fd.close()

parser = argparse.ArgumentParser(description='Extract amplicon performance stats from name sorted bam')
parser.add_argument('amplicon_bed', help='bed file of amplicon positions')
#parser.add_argument('reference', help='fasta reference in same coordinate system as amplicon bed')
parser.add_argument('bam', help='name sorted bam file input')
parser.add_argument('--filter', help='write out new bam file with annotated reads')

args = parser.parse_args()

def main():
    amplicon_set = args.amplicon_bed.split('.bed')[0].split('/')[-1]

    bam = pysam.AlignmentFile(args.bam, 'rb')

    filtered = None
    if args.filter:
        filtered = pysam.AlignmentFile(args.filter, 'wb', template=bam)

    amplicons = load_amplicons(args.amplicon_bed)
    histogram = bin_amplicons(amplicons, bam, filtered=filtered)
    total = sum(histogram.values())

    summary = '\t'.join([args.bam, str(total), args.amplicon_bed, args.reference])
    write_counts(sys.stdout, histogram)
    print(summary, sys.stderr)

if __name__ == "__main__":
    main()
