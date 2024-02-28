#!/bin/bash

############################################################################
##   Generate .cool and .mcool contact matrices from pairs using cooler   ##
##              https://cooler.readthedocs.io/en/latest/                  ##
############################################################################

### Activate conda environment before starting the script
#conda activate cooler_env 

### Supply external input file. 
### If you decide to use only interval user-defined variables, remove a piece of code below (if...fi)).

if [ -z "$6" ];
then echo -e "
You need to specify:
	1) path to output folder
	2) path to chromosome size file
	3) path to pairs files
	4) prefix output name of the .hic matrix
    5) number of CPU used
    6) bin size of the contact matrix in bp
	"
	exit
fi


### Input external and/or internal variables
OUTDIR="$1" 	    # specify your output directory (full path)
CHROMSIZE="$2"      # path to chromosome size file (full path)
INPUT="$3"  	    # path to the genome file in a fasta format  (full path)
f="$4"			    # prefix name of the output file
CPU="$5"            # number of CPU used
RES="$6"            # bin size
LOG="${OUTDIR}"/"${f}"_pairs2cool.log

### Check the arguments
echo -e \
"
Input arguments are
==========================
Output folder........: $1
Chromsize............: $2
Pairs................: $3
Prefix sample name...: $4
Number of CPU........: $5
Bin size in bp.......: $6
==========================

Make sure to activate conda environment before running this script!

" |& tee -a "${LOG}" 

### Make directories
mkdir -p "${OUTDIR}"/06_pairs2cool


### Load pairs into cooler format
echo -e \
"
===================================================
Load pairs into cool format with bin size ${RES} bp
===================================================
" |& tee -a "${LOG}"

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 "${CHROMSIZE}":"${RES}" \
    "${INPUT}" \
    "${OUTDIR}"/06_pairs2cool/"${f}"_"${RES}"bp.cool


### Balance the contact matrix using ICE
echo -e \
"
=============================================
Balance contact matrix with ICE (--mad-max 5)
=============================================
" |& tee -a "${LOG}"

cooler balance -p "${CPU}" --mad-max 5 --force "${OUTDIR}"/06_pairs2cool/"${f}"_"${RES}"bp.cool


### Generate multiresolution mcool files
echo -e \
"
====================================
Generate multiresolution mcool file
====================================
" |& tee -a "${LOG}"


cooler zoomify --nproc "${CPU}" \
               --balance \
               -o "${OUTDIR}"/06_pairs2cool/"${f}"_"${RES}"bp.mcool \
               "${OUTDIR}"/06_pairs2cool/"${f}"_"${RES}"bp.cool |& tee -a mlei_zoomify_CustomRes.log