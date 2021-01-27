"""Self-consistent alignment QC.
This is a wrapper around bamtools' pileup.
"""
import sys
import pysam

ref = pysam.FastaFile(sys.argv[2])
bam = pysam.AlignmentFile(sys.argv[1], 'rb')

for pc in bam.pileup(stepper='samtools', fastafile=ref):
    dels = 0
    ins = 0
    starts = 0
    ends = 0

    d = {'*': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0, '+': 0, '-': 0, ',': 0, '.': 0, 'N': 0}
    d_total = 0

    ref_base = ref.fetch(pc.reference_name, pc.reference_pos, pc.reference_pos+1)
    for read in pc.pileups:
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
            d[base] += 1
            d_total += 1

    freq = 0.0
    num_aligned = pc.get_num_aligned()
    ref_count = d[ref_base]
    if d_total == 0:
        #print(pc.reference_pos, "is empty")
        pass
    else:
        freq = ref_count / num_aligned

        if freq >= 0.95:
            continue

    bases = []
    for k in d:
        bases.append(f"{k}:{d[k]}")
    bases = ';'.join(bases)
    print('\t'.join(map(str, [pc.reference_name, pc.reference_pos + 1, ref_base, ref_count, d_total, pc.nsegments, num_aligned, freq, num_aligned / pc.nsegments, ins, dels, bases])))
