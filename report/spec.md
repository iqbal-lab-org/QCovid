## Summary of a dataset


A summary table:
95% of the samples have <N positions failing 85% criterion 
95% of the samples have <M positions with failing depth criterion
Y% of samples have no amplicons drops
Z% of samples have <=3 amplicon drops
ZZZ samples are completely excluded due to low quality


A single table, which is the concatenate of the per-sample tables. Columns are

1. overall depth
2. number of dropped amplicons (median depth across amplicon <5
3. Without excluding dropped amps, % of bases with <85% support
4. Without excluding dropped amps, % of bases with more than 5x depth and <85% support
5. Without excluding dropped amps, % of bases with <50% support
6. Without excluding dropped amps, % of bases with more than 5x depth and <50% support




In a future version we will add

7. After excluding dropped amps, % of bases with <85% support
8. After excluding dropped amps, % of bases with more than 5x depth and  <85% support
9. After excluding dropped amps, % of bases with <50% support
10. After excluding dropped amps, % of bases with more than 5x depth and <50% support

Two plots

1. y axis is count of samples. x axis is % of genome with >85% support and 5x depth (without excluding failed amplicons)
2. y axis is count of samples, x axis is dropped amplicons
3. histogram of amplicon numbers and counts of how often dropped.

## Summary of a technology (all nanopore, or all illumina)

As above. 

