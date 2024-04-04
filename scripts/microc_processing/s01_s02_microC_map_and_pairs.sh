#!/bin/bash

###########################################################
# Pipeline to preprocess and map to the genome Micro-C data
###########################################################

### Activate conda environment before starting the script
#conda activate cooler_env 

### Supply external input file. 
### If you decide to use only interval user-defined variables, remove a piece of code below (if...fi)).

if [ -z "$6" ];
then echo -e "
You need to specify:
	1) path to output folder
	2) path to genome file (fasta format)
	3) path to reads (1st pairs)
	4) path to reads (2nd pairs)
	5) number of threads
	6) prefix output name of the processed sample)
	"
	exit
fi

### Input external and internal variables
OUTDIR="$1" 	# specify your output directory (full path)
GENOME="$2" 	# path to the genome file in a fasta format  (full path)
FASTQ1="$3"		# read pair 1  (full path)
FASTQ2="$4"		# read pair 2  (full path)
CPU="$5"		# number of CPU used
f="$6"			# prefix name of the output file
MIN_RES=1000
LOG="${OUTDIR}"/"${f}".log

### Check the arguments
echo -e \
"
Input arguments are
==========================
Output folder......: $1
Genome.............: $2
Read 1.............: $3
Read 2.............: $4
Number of CPU......: $5
Prefix sample name.: $6
==========================

Make sure to activate conda environment before running this script!

" |& tee -a "${LOG}" 

### Make directories
mkdir -p "${OUTDIR}"/genome
mkdir -p "${OUTDIR}"/01_fastqc_report
mkdir -p "${OUTDIR}"/02_cooler_pipeline/bam
mkdir -p "${OUTDIR}"/02_cooler_pipeline/bw
mkdir -p "${OUTDIR}"/02_cooler_pipeline/pairs
mkdir -p "${OUTDIR}"/02_cooler_pipeline/stats
mkdir -p "${OUTDIR}"/02_cooler_pipeline/cooler/cload
mkdir -p "${OUTDIR}"/02_cooler_pipeline/cooler/mcool

### Prepare bwa genome index 
echo -e \
"Prepare genome index with bwa
=============================="  |& tee -a "${LOG}"

cd "${OUTDIR}"/genome || exit
base="$(basename "${GENOME}" )"
ln -s "${GENOME}" "${base}"

bwa index "${base}"  |& tee -a "${LOG}"
INDEX="${OUTDIR}"/genome/"${base}"

### Prepare chromosome size file
samtools faidx "${base}"
cut -f1,2 "${base}".fai > "${base}".chromsize
CHROMSIZE="${OUTDIR}"/genome/"${base}".chromsize

### Prepare genomic bins of fixed size
cooler makebins "${CHROMSIZE}" 200 > "${base}".200bp_cooler.bed
DIGEST200="${OUTDIR}"/genome/"${base}".200bp_cooler.bed

### FastQC quality check of the input reads
echo -e \
"
==========================================
FastQC quality check of the input reads
==========================================
" |& tee -a "${LOG}"

cd "${OUTDIR}"/01_fastqc_report || exit
fastqc -o ./ "${FASTQ1}" -t "${CPU}" |& tee -a "${LOG}"
fastqc -o ./ "${FASTQ2}" -t "${CPU}" |& tee -a "${LOG}"

### Map raw reads and process them into pairs and cool files
echo -e \
"
==========================================
Map raw reads and process them into pairs
==========================================
" |& tee -a "${LOG}"

bwa mem -SP5M -t"${CPU}" \
	"${INDEX}" \
	"${FASTQ1}" \
	"${FASTQ2}" \
| pairtools parse \
	--chroms-path "${CHROMSIZE}" \
	--walks-policy 5unique \
	--max-inter-align-gap 20 \
    --nproc-in "${CPU}" --nproc-out "${CPU}" \
| pairtools sort  \
	--nproc "${CPU}" \
| pairtools dedup \
	--nproc-in "${CPU}" --nproc-out "${CPU}" \
	--mark-dups  \
	--max-mismatch 2 \
 	--output-stats "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt \
	--output-unmapped "${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_unmapped.pairs.gz \
	--output-dups "${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_dups.pairs.gz \
| pairtools split \
	--nproc-in "${CPU}" --nproc-out "${CPU}" \
	--output-pairs "${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_nodups.pairs.gz  \
	--output-sam - \
| samtools view -bS - \
| samtools sort -@"${CPU}" - -o "${OUTDIR}"/02_cooler_pipeline/bam/"${f}"_nodups.bam

samtools index "${OUTDIR}"/02_cooler_pipeline/bam/"${f}"_nodups.bam

bamCoverage -b "${OUTDIR}"/02_cooler_pipeline/bam/"${f}"_nodups.bam \
	-bs 10 -p "${CPU}" --ignoreDuplicates \
	-o "${OUTDIR}"/02_cooler_pipeline/bw/"${f}"_nodups.bw

# Prepare the txt input to estimate the percentage of cis and trans contacts in R
	total=$(grep 'total_nodups' "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt | awk '{print $0}')
	cis=$(grep 'cis\>' "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt | awk '{print $0}')
	trans=$(grep 'trans' "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt | awk '{print $0}')
	cis_1kb=$(grep 'cis_1kb+' "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt | awk '{print $0}')
	cis_10kb=$(grep 'cis_10kb+' "${OUTDIR}"/02_cooler_pipeline/stats/"${f}"_pairtools.stats.txt | awk '{print $0}')
	{
		echo -e "${total}"
		echo -e "${cis}"
		echo -e "${trans}"
		echo -e "${cis_1kb}"	
		echo -e "${cis_10kb}"
	} >> "${OUTDIR}"/02_cooler_pipeline/stats/"${f}".distance.txt


# Filter self-ligated pairs or pairs from close nucleosomes (= remove reads that fall into the same 200 bp bins).
pairtools restrict -f "${DIGEST200}" \
	--nproc-in "${CPU}" --nproc-out "${CPU}" \
	"${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_nodups.pairs.gz \
| pairtools select \
	--nproc-in "${CPU}" --nproc-out "${CPU}" \
	-o "${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_nodups.same200bp_bin.pairs.gz \
	--output-rest "${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_nodups.200bp_bins_filtered.pairs.gz \
	'(chrom1==chrom2) and (COLS[9]==COLS[12])'

# Load pairs into cooler format at a resolution MIN_RES (1000 bp)
echo -e \
"
=============================
Load pairs into cooler format
==============================
" |& tee -a "${LOG}"

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "${CHROMSIZE}":${MIN_RES} \
	"${OUTDIR}"/02_cooler_pipeline/pairs/"${f}"_nodups.200bp_bins_filtered.pairs.gz \
	"${OUTDIR}"/02_cooler_pipeline/cooler/cload/"${f}"_${MIN_RES}Cload.cool |& tee -a "${LOG}"


# Balance cool matrices
echo -e \
"
=============================================
Balance contact matrix with ICE (--mad-max 5)
=============================================
" |& tee -a "${LOG}"

cooler balance -p "${CPU}" --mad-max 5 --force \
	"${OUTDIR}"/02_cooler_pipeline/cooler/cload/"${f}"_${MIN_RES}Cload.cool |& tee -a "${LOG}"


# Generate multiresolution mcool files
echo -e \
"
====================================
Generate multiresolution mcool file
====================================
" |& tee -a "${LOG}"

cooler zoomify --nproc "$CPU" --balance -r 1000,2000,5000,10000,25000,50000,100000 \
	-o "${OUTDIR}"/02_cooler_pipeline/cooler/mcool/"${f}"_${MIN_RES}Cload.mcool \
	"${OUTDIR}"/02_cooler_pipeline/cooler/cload/"${f}"_${MIN_RES}Cload.cool |& tee -a "${LOG}"

echo -e \
"
====================================
                Done!
====================================
" |& tee -a "${LOG}"