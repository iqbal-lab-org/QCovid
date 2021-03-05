"""bin reads by amplicon coordinates
"""

import sys
import pysam
import numpy as np

# loosen the coordinate interval to cast non-softclipped adapter sequences
PADDING = 5

global_stats = {'reads': 0, 'unmapped': 0, 'unmated': 0, 'secondary': 0, 'mate_off_target': 0, 'overlap_primer_boundary': 0, 'r1_primer_match_only': 0, 'fragment_size_mean': 0.0, 'fragment_size_std': 0.0, 'read_length_mean': 0.0, 'read_length_std': 0.0, 'F1R2': 0, 'R1F2': 0, 'F1F2': 0, 'R1R2': 0}

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
    total = 0
    for amplicon in amplicons:
        start, end = amplicons[amplicon]

        # ensure coords are sorted
        rc = 0
        for read in bam.fetch(reference, start, end):
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

def bam_stats(bam):
    gs = {'reads': 0, 'unmapped': 0, 'unmated': 0, 'secondary': 0, 'mate_off_target': 0, 'overlap_primer_boundary': 0, 'r1_primer_match_only': 0, 'fragment_size_mean': 0.0, 'fragment_size_std': 0.0, 'read_length_mean': 0.0, 'read_length_std': 0.0, 'r1forward': 0, 'r1reverse': 0, 'r2forward': 0, 'r2reverse': 0}
    tlens = []
    rlens = []

    for aln in bam.fetch():

        gs['reads'] += 1

        if aln.is_unmapped:
            gs['unmapped'] += 1

        if aln.mate_is_unmapped:
            gs['unmated'] += 1

        if aln.is_secondary:
            gs['secondary'] += 1
    
        if aln.is_read1 and aln.is_reverse:
            gs['r1reverse'] += 1
        if aln.is_read1 and not aln.is_reverse:
            gs['r1forward'] += 1
        if aln.is_read2 and not aln.is_reverse:
            gs['r2forward'] += 1
        if aln.is_read2 and aln.is_reverse:
            gs['r2reverse'] += 1

        #if aln.is_read1: # and not aln.is_reverse:
        tlens.append(aln.template_length)
        rlens.append(aln.query_length)

    gs['fragment_size_mean'] = np.mean(np.abs(tlens))
    gs['fragment_size_std'] = np.std(np.abs(tlens))

    gs['read_length_mean'] = np.mean(np.abs(rlens))
    gs['read_length_std'] = np.std(np.abs(rlens))

    print(max(tlens))
    print(min(tlens))
    return gs


def write_counts(amps_fd, histogram):
    for a in histogram:
        print('\t'.join([a, str(histogram[a])]), file=amps_fd)
    amps_fd.close()

def main():
    amplicon_bed = sys.argv[1]
    reference = "MN908947"

    amplicon_set = amplicon_bed.split('.bed')[0].split('/')[-1]

    basename = sys.argv[2].split('.')[0]

    out_name = f"{basename}.{amplicon_set}.bam"
    out_amps = f"{basename}.{amplicon_set}.tsv"
    out_stats = f"{basename}.stats.tsv"

    bam = pysam.AlignmentFile(sys.argv[2], 'rb')
    #filtered = pysam.AlignmentFile(out_name, 'wb', template=bam)
    filtered = None
    amplicons = load_amplicons(amplicon_bed)
    histogram = bin_amplicons(reference, amplicons, bam, filtered=filtered)
    total = sum(histogram.values())

    amps_fd = open(out_amps, 'w')
    summary = '\t'.join([basename, str(total), amplicon_bed, reference])
    write_counts(amps_fd, histogram)
    amps_fd.close()
    bam.close()
    print(summary)

    bam = pysam.AlignmentFile(sys.argv[2], 'rb')
    bs = bam_stats(bam)
    bam.close()

    outf = open(out_stats, 'w')
    print('\t'.join(list(map(str, bs.keys()))), file=outf)
    print('\t'.join(list(map(str, bs.values()))), file=outf)
    outf.close()

if __name__ == "__main__":
    main()
