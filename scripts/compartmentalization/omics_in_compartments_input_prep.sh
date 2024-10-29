#!/bin/bash

# For each ATAC-seq and ChIP-seq experiment, calculate the CPM-normalized coverage using a window 
# size that results in 20,000 bins across the genome. Below is an example for the H3K4me3 ChIP-seq signal in M. leidyi.

bamCoverage -b ./Mlei_H3K4me3.bam \
	--numberOfProcessors 16 \
	-o mlei_chip_H3K4me3_20Kbins.bed -of bedgraph \
	--binSize 10252 \
	--normalizeUsing CPM \
	--extendReads \
	--ignoreDuplicates

