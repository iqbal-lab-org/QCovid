import sys
import pysam

amplicons = {}
amplicon_bed = "nCoV-nl-primal500-75.bed"
reference = "MN908947.3"

for line in open(amplicon_bed):
    if line[0] == '#':
        continue
    fname, start, end = line.strip().split()
    amplicons[fname] = [int(start), int(end)]

#print(amplicons)
bam = pysam.AlignmentFile(sys.argv[1], 'rb')
basename = sys.argv[1].split('.')[0]
amplicounts = {}

out_name = f"{basename}.amplicons.bam"
out_amps = f"{basename}.amplicons.tsv"

filtered = pysam.AlignmentFile(out_name, 'wb', template=bam)
total = 0
for amplicon in amplicons:
    start, end = amplicons[amplicon]

    # ensure coords are sorted
    rc = 0
    for read in bam.fetch(reference, start, end):
        if read.reference_start > (start - 5) and read.reference_end < (end + 5):
            rc += 1 
            filtered.write(read)
            total += 1
            if amplicon == 'SARS-CoV-2_36_pool2':
                print(read.reference_start, read.reference_end, read.query_alignment_start, read.query_length, read.query_alignment_sequence)
    amplicounts[amplicon] = rc

filtered.close()
bam.close()

amps_fd = open(out_amps, 'w')
print('\t'.join([basename, str(total), amplicon_bed, reference]))
for a in amplicounts:
    print('\t'.join([a, str(amplicounts[a])]), file=amps_fd)
amps_fd.close()
