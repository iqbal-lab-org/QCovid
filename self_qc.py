"""Self-consistent alignment QC.
This is a wrapper around bamtools' pileup.
"""
import sys
import pysam

ref = pysam.FastaFile(sys.argv[2])
bam = pysam.AlignmentFile(sys.argv[1], 'rb')
fastout = '--fasta' in sys.argv

if not fastout:
    print('\t'.join(["reference", "position", "ref", "ref_depth", "total_depth", "n_segments", "n_aligned", "freq", "n_segments/n_aligned", "insertions", "deletions", "bases"]))

refseq = ""
orig_name = ""
for line in open(sys.argv[2]):
    if line[0] == '>':
        orig_name = line[1:]
        continue
    refseq += line.strip()
masked = list(refseq)

for pc in bam.pileup(stepper='samtools', fastafile=ref):
    dels = 0
    ins = 0
    starts = 0
    ends = 0

    d = {'*': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0, '+': 0, '-': 0, ',': 0, '.': 0, 'N': 0}
    d_total = 0

    ref_base = ref.fetch(pc.reference_name, pc.reference_pos, pc.reference_pos+1)
    if ref_base not in d:
        ref_base = 'N'
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
            if base not in d:
                base = 'N'
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
    if not fastout:
        print('\t'.join(map(str, [pc.reference_name, pc.reference_pos + 1, ref_base, ref_count, d_total, pc.nsegments, num_aligned, freq, num_aligned / pc.nsegments, ins, dels, bases])))
    else:
        masked[pc.reference_pos] = 'N'

if fastout:
    print(f">{orig_name}_masked")
    print(''.join(masked))
