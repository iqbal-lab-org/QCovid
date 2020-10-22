bwa mem MN908947.fasta $3 | samtools view -bS - > $WORKDIR/$2.MN908947.bam
samtools sort $WORKDIR/$2.MN908947.bam -o $WORKDIR/$2.MN908947.sorted.bam
samtools index $WORKDIR/$2.MN908947.sorted.bam
python3 bin_amplicons.py primers/$1.bed $WORKDIR/$2.MN908947.sorted.bam
samtools sort $WORKDIR/$2.$1.bam -o $WORKDIR/$2.$1.sorted.bam
samtools index $WORKDIR/$2.$1.sorted.bam
