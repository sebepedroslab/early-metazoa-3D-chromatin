# **Micro-C processing pipeline**

1. Align raw reads to the reference genome, filter valid alignments, parse mapped reads into pairs files followed by matrix aggregation and normalization.
>
>[s01_s02_microC_map_and_pairs.sh](s01_s02_microC_map_and_pairs.sh)
>

2. Estimate reproducibility between replicates using the stratum-adjusted correlation coefficient (SCC).
>
>[s03_hicRep_replicate_correlation.sh](s03_hicRep_replicate_correlation.sh)
>

3. Merge correlated biological and technical replicates into a single pairs file.
>
>[s04_mergeReplicates.sh](s04_mergeReplicates.sh)
>

4. Aggregate valid pairs into normalized multiresolution contact matrix.
>Contact matrices were normalized with either juicer tools using using Knight-Ruiz balancing method
>[s05_pairs2hic.sh](s05_pairs2hic.sh)<br>
or with cooler using ICE balancing method<br>
>[s06_pairs2cool.sh](s06_pairs2cool.sh)

