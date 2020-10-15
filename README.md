# QCovid
QC pipelines for sars-cov-2 sequence+consensus submitted to the ENA

This QC pipeline considers raw sequencing data and is flexible over:

* Single vs. paired reads (i.e. Nanopore vs. Illumina)
* Amplicon primers

## Single and Paired read pipelines

We use `bwa` to map reads as part of the pipeline. Downstream analysis of the mapping is agnostic to mate pairing.

## Preprocessing primer sets

Primers should be presented as a list of tab delimited names, start positions, and end positions.

e.g. a bed file:

```
#Reference: MN908947.fasta (29903 bp.)
SARS-CoV-2_1_pool1	30	495
SARS-CoV-2_3_pool1	704	1205
SARS-CoV-2_5_pool1	1372	1897
SARS-CoV-2_7_pool1	2115	2642
SARS-CoV-2_9_pool1	2786	3257
SARS-CoV-2_11_pool1	3460	3948
SARS-CoV-2_13_pool1	4111	4634
SARS-CoV-2_15_pool1	4809	5359
SARS-CoV-2_17_pool1	5563	6037
```

There are about a dozen primer sets. Most samples post-April are ARTIC-v3.

Primer sets are not defined consistently and have to be preprocessed into these bed files individually. In general this entails parsing some list of sequences and putting them into a consistent coordinate system. We are using MN908947. If possible the amplicon names should follow the convention `SARS-CoV-2_{number}_{pool}` where `number` is the ID of the amplicon and `pool` specifies which (if any) reaction mix the primer is included in.

For reference these scripts and raw primer sets are included in `primers/`. Brief notes follow.

### Dutch Primal 500-75
`nCoV-nl-primal500-75.bed`. These were extracted from the supplementary materials for [...], provided in `primers/dutch_primal_amplicons.txt` and preprocessed with the script `primers/parse_nl_amplicons.py`.

### ARTIC-v3
`nCoV-artic-v3.bed`.

## Fetching raw sequencing data

`ena_data_get`

## Target datasets

### Nanopore

* PRJEB37966
* PRJEB37886  https://www.cogconsortium.uk/protocols/
* PRJEB38388  https://www.biorxiv.org/content/10.1101/2020.04.21.050633v1.supplementary-material
* PRJNA613958 https://peerj.com/articles/9255/
* PRJNA614995 https://peerj.com/articles/9255/
* PRJNA616147 ARTIC v1 primers
* PRJNA622817 https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html
* PRJNA627229 https://www.medrxiv.org/content/10.1101/2020.04.25.20079517v1.full.pdf
* PRJNA610248
* PRJNA614976 https://doi.org/10.1016/j.cell.2020.04.021
* PRJNA632678 https://doi.org/10.1128/MRA.00573-20
* PRJNA601630 https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext
* PRJNA608224 https://www.biorxiv.org/content/10.1101/2020.03.05.976167v1.full
* PRJNA608242
* PRJNA612578 https://www.biorxiv.org/content/biorxiv/early/2020/04/12/2020.04.09.034462.full.pdf
* PRJNA614504 https://europepmc.org/article/med/32676620
* PRJNA627286 https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.26140

### Illumina
* TODO

## Running the pipeline

For single-end reads and the Dutch primers (PRJxxxxxxxxx), one would run:

```bash
bwa mem MN908947.fasta $PROJECT/$SAMPLE_1.fastq.gz | samtools view -bS - > $WORK/$SAMPLE.MN908947.bam
samtools sort $WORK/$SAMPLE.MN908947.bam -o $WORK/$SAMPLE.MN908947.sorted.bam
samtools index $WORK/$SAMPLE.MN908947.sorted.bam
python3 bin_amplicons.py primers/nCoV-nl-primal500-75.bed $SAMPLE.MN908947.sorted.bam
samtools sort $WORK/$SAMPLE.amplicons.bam -o $WORK/$SAMPLE.amplicons.sorted.bam
samtools index $WORK/$SAMPLE.amplicons.sorted.bam
````

TODO: include run-time estimates, nextflow batching by project

## Interpreting results

### Collating amplicon depths across a project

TODO: jupyter notebook heatmap

### Per-sample internal consistency

### Investigating cleaned reads
`bin_amplicons.py` will optionally write reads that positively identify as amplicons to a 'cleaned' sam file.

