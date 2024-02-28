#!/bin/bash

#####################################################################################
# To merge biological and technical replicates into one pairs files, we used a script
# from the Docker-4dn-hic repository:
# https://github.com/4dn-dcic/docker-4dn-hic/blob/master/scripts/run-merge-pairs.sh
#####################################################################################

### Activate conda environment before starting the script
#conda activate cooler_env 


### Make directories
mkdir -p /media/data/04_merged_pairs


### Define the following variables
CPU=10													# specify number of CPUs (only important if you're going to merge technical replicates)
INDIR="/media/data/02_cooler_pipeline/pairs"			# specify your input directory with pairs files (full path)
OUTDIR="/media/data/04_merged_pairs"					# specify your output directory (full path)
SCRIPT_4DN=/media/data/docker-4dn-hic-master/scripts	# path to the directory with 4DN scripts
MERGED_PAIRS="sample"									# file name of merged pairs


### Check the arguments
echo -e \
"
Input arguments are
=============================
Output folder......: $OUTDIR
Input pairs folder.: $INDIR
=============================

Make sure to activate conda environment before running this script!
"

cd ${OUTDIR} || exit

### To merge technical replicates, you first need to remove duplicated reads. 
### These are PCR duplicates produced from the same duplication event.

TechRep1=${INDIR}/"technical_replicate_1.pairs.gz"		# specify the path to the technical replicate 1
TechRep2=${INDIR}/"technical_replicate_2.pairs.gz"		# specify the path to the technical replicate 2

pairtools merge \
	${TechRep1} \
	${TechRep2} \
| pairtools dedup \
	--nproc-in ${CPU} --nproc-out ${CPU} \
	--mark-dups  \
	--max-mismatch 2 \
 	--output-stats ${OUTDIR}/technical_merged_dedup_pairtools.stats.txt \
	--output-dups ${OUTDIR}/technical_merged_duplicated.pairs.gz \
	--output ${OUTDIR}/technical_merged.pairs.gz


### Merge biological replicates and deduplicated technical replicates that pass 
### the SCC correlation analysis (f03_HicRep_replicate_correlation.sh) using 4DN script:

BioRep1=${INDIR}/"biological_replicate_1.pairs.gz"		# specify the path to the biological replicates to be merged
BioRep2=${INDIR}/"biological_replicate_2.pairs.gz"		# specify the path to the biological replicates to be merged

${SCRIPT_4DN}/run-merge-pairs.sh ${MERGED_PAIRS} \
  ${BioRep1} \
  ${BioRep2} \
  ${OUTDIR}/technical_merged.pairs.gz


### Generate stats file of merged pairs
pairtools stats \
	-o ${OUTDIR}/${MERGED_PAIRS}.stats \
	${OUTDIR}/${MERGED_PAIRS}.pairs.gz