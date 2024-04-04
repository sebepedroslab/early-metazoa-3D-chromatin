#!/bin/bash

############################################################################
##   Generate .hic format contact matrix from pairs using juicer tools    ##
##              https://github.com/aidenlab/juicer/wiki/Pre               ##
############################################################################

### Supply external input file. 
### If you decide to use only interval user-defined variables, remove a piece of code below (if...fi)).

if [ -z "$5" ];
then echo -e "
You need to specify:
	1) path to output folder
	2) path to chromosome size file
	3) path to pairs files
	4) prefix output name of the .hic matrix
    5) path to the juicer_tools.jar 
	"
	exit
fi


### Input external and/or internal variables
OUTDIR="$1" 	    # specify your output directory (full path)
CHROMSIZE="$2"      # path to chromosome size file (full path)
INPUT="$3"  	    # path to the genome file in a fasta format  (full path)
f="$4"			    # prefix name of the output file
JUICER_TOOLS="$5"   # path to the juicer_tools.jar
LOG="${OUTDIR}"/"${f}"_pairs2hic.log


### Check the arguments
echo -e \
"
Input arguments are
==========================
Output folder........: $1
Chromsize............: $2
Pairs................: $3
Prefix sample name...: $4
Path to juicer_tools.: $5
==========================

" |& tee -a "${LOG}" 

### Make directories
mkdir -p "${OUTDIR}"/05_pairs2hic


### Generate contact matrix .hic 
java -Xmx200g -jar "${JUICER_TOOLS}"/juicer_tools.jar pre \
        "${INPUT}" "${OUTDIR}"/05_pairs2hic/"${f}".hic \
        "${CHROMSIZE}" \
        -k VC,VC_SQRT,KR,SCALE \
        -r 200,400,800,1000,2000,4000,5000,7000,10000,25000,50000,100000


### Validate .hic file
java -Xmx200g -jar "${JUICER_TOOLS}"/juicer_tools.jar validate "${OUTDIR}"/05_pairs2hic/"${f}".hic