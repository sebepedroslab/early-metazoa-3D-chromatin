#!/bin/bash

###################################################################################
# To identify the resolution and window sizes reflecting the strongest partitioning
# of the target genome into isolated domains, we generated aggregated metaplots of
# insulation boundaries. The resolution and two window sizes with maximum average 
# insulation scores were consirered to be optimal and used for further analysis.
###################################################################################

### Define the following variables:
OUTDIR="./data/insulation"
CHROMSIZE="Mlei_gDNA.chromsize"
species="mlei"

### Create output directories
mkdir -p ${OUTDIR}/${species}_bg
mkdir -p ${OUTDIR}/metaplot

### Prepare input data for metaplots
# Select genomic coordinates of insulation boundaries marked as strong using Li threshold.
for f in ${species}_400bpRes ${species}_1000bpRes ${species}_2000bpRes ${species}_4000bpRes; do
awk -v OFS='\t' '($5=="False" && ($16=="True" || $17=="True" || $18=="True")) {print $1"\t"$2"\t"$3}' ${OUTDIR}/"${f}".boundaries_Li.tsv > \
${OUTDIR}/"${f}"_strong_boundary_region.bed

# Extend the region around insulation boundary to 15 kb upstream and downstream using bedtools.
bedtools slop -i ${OUTDIR}/"${f}"_strong_boundary_region.bed \
              -g ${CHROMSIZE} \
              -b 15000 > \
              ${OUTDIR}/metaplot/"${f}"_strong_boundary_region_slop.bed

# Split each extended interval into 30 windows and label windows by window number using bedtools.
bedtools makewindows -b ${OUTDIR}/metaplot/"${f}"_strong_boundary_region_slop.bed \
                     -n 60 -i winnum  > \
                     ${OUTDIR}/metaplot/"${f}"_strong_boundary_region_slop_windows.bed

rm ${OUTDIR}/metaplot/"${f}"_strong_boundary_region_slop.bed

# Generate BedGraph files with insulation score for different sliding windows and resolutions
awk -v OFS='\t' '($5=="False") {print $1"\t"$2"\t"$3"\t"$7"\t""wx5"}' ${OUTDIR}/"${f}".boundaries_Li.tsv > ${OUTDIR}/${species}_bg/"${f}"_w+s_region.bg
awk -v OFS='\t' '($5=="False") {print $1"\t"$2"\t"$3"\t"$9"\t""wx10"}' ${OUTDIR}/"${f}".boundaries_Li.tsv >> ${OUTDIR}/${species}_bg/"${f}"_w+s_region.bg
awk -v OFS='\t' '($5=="False") {print $1"\t"$2"\t"$3"\t"$11"\t""wx25"}' ${OUTDIR}/"${f}".boundaries_Li.tsv >> ${OUTDIR}/${species}_bg/"${f}"_w+s_region.bg
done