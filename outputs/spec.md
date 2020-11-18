
1. reports, described elsewhere
2. single index file (ideally SQL database) which lists all samples and includes an entry for include or exclude and a reason.
3. for each sample, masked consensus for high quality (dropped amp and positions with <85% support, and depth< (10 for illumina and 10 for nanopore))
4. for each sample, masked consensus excluding dropped amplicons.
5. for each sample, 4 masks as bed files (85%, amplicons, depth)
6. unmasked input consensus
