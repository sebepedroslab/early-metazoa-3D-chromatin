#!/bin/bash

#################################################################################
# SCC correlation across replicates at different map resolutions (python version)
# doi: 10.1093/bioinformatics/btab097; doi: 10.1101/gr.220640.117
#################################################################################

### Activate conda environment before starting the script
#conda activate cooler_env 

### Put all mcool files of replicates being compared in the same directory


### Define the following variables
OUTDIR="/media/data/test" 	                                    # specify your output directory (full path)
INPUT="/media/data/test/02_cooler_pipeline/cooler/mcool" 	      # path to the mcool files (full path)
LOG="${OUTDIR}"/f03_SCC_replicate_correlation.log                 # path to the log output file

### Check the arguments
echo -e \
"
Input arguments are
=============================
Output folder......: $OUTDIR
Input mcool folder.: $INPUT
=============================

Make sure to activate conda environment before running this script!

" |& tee -a "${LOG}" 

### Make directories
mkdir -p "${OUTDIR}"/03_compare_replicates/res1000
mkdir -p "${OUTDIR}"/03_compare_replicates/res2000
mkdir -p "${OUTDIR}"/03_compare_replicates/res5000
mkdir -p "${OUTDIR}"/03_compare_replicates/res10000
mkdir -p "${OUTDIR}"/03_compare_replicates/res25000
mkdir -p "${OUTDIR}"/03_compare_replicates/res50000

cd "${INPUT}" || exit

### calculate for resolution 1000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 1,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
       hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res1000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 1000 --h 10 --dBPMax 5000000
	 hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res1000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 1000 --h 10 --dBPMax 5000000
	 hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res1000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 1000 --h 10 --dBPMax 5000000;
     done;
done 


### calculate for resolution 2000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 2,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
      hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res2000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 2000 --h 10 --dBPMax 5000000
	hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res2000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 2000 --h 10 --dBPMax 5000000
	hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res2000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 2000 --h 10 --dBPMax 5000000;
     done;
done 


### calculate for resolution 5000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 5,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
      hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res5000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 5000 --h 10 --dBPMax 5000000
	hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res5000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 5000 --h 10 --dBPMax 5000000
	hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res5000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 5000 --h 10 --dBPMax 5000000;
      done;
done 


### calculate for resolution 10000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 10,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
      hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res10000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 10000 --h 10 --dBPMax 5000000
	hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res10000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 10000 --h 10 --dBPMax 5000000
	hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res10000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 10000 --h 10 --dBPMax 5000000;
      done;
done 


### calculate for resolution 25000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 25,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
      hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res25000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 25000 --h 10 --dBPMax 5000000
	hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res25000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 25000 --h 10 --dBPMax 5000000
	hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res25000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 25000 --h 10 --dBPMax 5000000;
      done;
done 


### calculate for resolution 50000 no Downsampling
echo -e \
"
==========================================
        SCC correlation at 50,000 bp
==========================================
" |& tee -a "${LOG}"

f=(*)
for ((i = 0; i < ${#f[@]}; i++)); do 
      for ((j = i + 1; j < ${#f[@]}; j++)); do 
      hicrep "${f[i]}" "${f[j]}" "${OUTDIR}"/03_compare_replicates/res50000/"${f[i]}"+"${f[j]}"_output.SCC --binSize 50000 --h 10 --dBPMax 5000000
	hicrep "${f[i]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res50000/"${f[i]}"+"${f[i]}"_output.SCC --binSize 50000 --h 10 --dBPMax 5000000
	hicrep "${f[j]}" "${f[i]}" "${OUTDIR}"/03_compare_replicates/res50000/"${f[j]}"+"${f[i]}"_output.SCC --binSize 50000 --h 10 --dBPMax 5000000;
      done;
done 


echo -e \
"
==========================================
                   Done!
==========================================
" |& tee -a "${LOG}"