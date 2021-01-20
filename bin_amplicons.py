"""bin reads by amplicon coordinates
"""

import sys
import pysam

def load_amplicons(amplicon_bed):
    amplicons = {}
    for line in open(amplicon_bed):
        if line[0] == '#':
            continue
        fname, start, end = line.strip().split()
        amplicons[fname] = [int(start), int(end)]
    return amplicons

def bin_amplicons(reference, amplicons, bam, filtered=None):
    histogram = {}
    for amplicon in amplicons:
        start, end = amplicons[amplicon]

        # ensure coords are sorted
        rc = 0
        for read in bam.fetch(reference, start, end):
            if not read.reference_start or not read.reference_end:
                continue
            if read.reference_start > (start - 5) and read.reference_end < (end + 5):
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

def main():
    amplicon_bed = sys.argv[1]
    reference = "MN908947.3"

    amplicon_set = amplicon_bed.split('.bed')[0].split('/')[-1]

    basename = sys.argv[2].split('.')[0]

    out_name = f"{basename}.{amplicon_set}.bam"
    out_amps = f"{basename}.{amplicon_set}.tsv"

    bam = pysam.AlignmentFile(sys.argv[2], 'rb')
    filtered = pysam.AlignmentFile(out_name, 'wb', template=bam)
    amplicons = load_amplicons(amplicon_bed)
    histogram = bin_amplicons(reference, amplicons, bam, filtered=filtered)
    total = sum(histogram.values)

    amps_fd = open(out_amps, 'w')
    summary = '\t'.join([basename, str(total), amplicon_bed, reference])
    write_counts(amps_fd, histogram)
    amps_fd.close()
    print(summary)

if __name__ == "__main__":
    main()
