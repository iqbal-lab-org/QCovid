WORKDIR=work
python3 qcovid.py $1 fasta $2 > $WORKDIR/$2.fasta
minimap2 -ax sr $WORKDIR/$2.fasta $3 $4 | samtools sort -o $WORKDIR/$2.bam
samtools index $WORKDIR/$2.bam
python3 self_qc.py $WORKDIR/$2.bam $WORKDIR/$2.fasta --fasta > $WORKDIR/$2.masked.fasta
python3 self_qc.py $WORKDIR/$2.bam $WORKDIR/$2.fasta > $WORKDIR/$2.self_qc.tsv
python3 qcovid.py $1 import --self_qc $WORKDIR/$WORKDIR/$2.self_qc.tsv
