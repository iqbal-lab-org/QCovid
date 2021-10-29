## Primer sets

### Dutch Primal 500-75
`nCoV-nl-primal500-75.bed`. These were extracted from the supplementary materials for [...], provided in `primers/dutch_primal_amplicons.txt` and preprocessed with the script `primers/parse_nl_amplicons.py`. Note that some of these primers were supplied with IUPAC ambiguous nucleotides which may present an issue with locating them in the reference sequence.

### ARTIC-v3

Source: https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv

`python parse_articv3_amplicons.py > nCoV-artic-v3.bed`

### ARTIC-v4

Source: https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V4/SARS-CoV-2.primer.bed

### Midnight-1200

Source: https://docs.google.com/spreadsheets/d/1M5I_C56ZC8_2Ycgm9EFieVlVNqxsP7dXAnGoBZy3nDo/edit#gid=755704891


## Target datasets

### Nanopore

* `PRJEB37966`
* `PRJEB37886` ARTIC v3(?) https://www.cogconsortium.uk/protocols/
* `PRJEB38388` Dutch primers https://www.biorxiv.org/content/10.1101/2020.04.21.050633v1.supplementary-material
* `PRJNA613958` https://peerj.com/articles/9255/
* `PRJNA614995` https://peerj.com/articles/9255/
* `PRJNA616147` ARTIC v1 primers
* `PRJNA622817` https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html
* `PRJNA627229` https://www.medrxiv.org/content/10.1101/2020.04.25.20079517v1.full.pdf
* `PRJNA610248`
* `PRJNA614976` https://doi.org/10.1016/j.cell.2020.04.021
* `PRJNA632678` https://doi.org/10.1128/MRA.00573-20
* `PRJNA601630` https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30154-9/fulltext
* `PRJNA608224` https://www.biorxiv.org/content/10.1101/2020.03.05.976167v1.full
* `PRJNA608242`
* `PRJNA612578` https://www.biorxiv.org/content/biorxiv/early/2020/04/12/2020.04.09.034462.full.pdf
* `PRJNA614504` https://europepmc.org/article/med/32676620
* `PRJNA627286` https://onlinelibrary.wiley.com/doi/full/10.1002/jmv.26140

### Illumina
* TODO
