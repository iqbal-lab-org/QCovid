bwa mem MN908947.fasta $2 $3 | samtools view -bS - > $WORKDIR/$1.MN908947.bam
samtools sort $WORKDIR/$1.MN908947.bam -o $WORKDIR/$1.MN908947.sorted.bam
samtools index $WORKDIR/$1.MN908947.sorted.bam
python3 bin_amplicons.py primers/nCoV-generic-300bp-windows.bed $WORKDIR/$1.MN908947.sorted.bam
samtools sort $WORKDIR/$1.nCoV-generic-300bp-windows.bam -o $WORKDIR/$1.nCoV-generic-300bp-windows.sorted.bam
samtools index $WORKDIR/$1.nCoV-generic-300bp-windows.sorted.bam

