#!/bin/bash

###############################################################################################
# For each species, we calculated insulation scores profiles at resolutions approximately 
# equivalent to 50,000, 100,000, 200,000 and 400,000 genomic bins per species genome with
# a sliding window for each resolution that is 5x, 10x and 25x the applied resolution.
# The applied resolutions corresponded to bin sizes of 400 bp, 1,000 bp, 2,000 bp and 4,000 bp
# for M. leidyi and E. muelleri; for S. arctica we applied resolutions of 400bp, 800bp, 1,600bp,
# and 2,800 bp; for C. owczarzaki, S. rosetta, and T. adhaerens, insulation scores were called
# at resolutions 200 bp, 400 bp, 800 bp, and 1,600 bp; for N. vectensis - 500 bp, 1,000 bp, 
# 2,500 bp, 5,000 bp. Below is an example of the script used to calculate insultaion profiles 
# in Mnemiopsis leidyi.
###############################################################################################

### Activate conda environment before starting the script
conda activate cooltools_env

### Make directories
mkdir -p ./data/insulation

### Define the following variables:
OUTDIR="./data/insulation"
mcool="mlei_200bp.mcool"
species="mlei"

### Calculate insulation profiles at different resolutions and window sizes
# Resolution 400 bp
cooltools insulation -o ${OUTDIR}/${species}_400bpRes.boundaries_Li.tsv \
                     --threshold 'Li' ${mcool}::/resolutions/400 2000 4000 10000  \
                     --bigwig --min-dist-bad-bin 2

# Resolution 1,000 bp
cooltools insulation -o ${OUTDIR}/${species}_1000bpRes.boundaries_Li.tsv \
                     --threshold 'Li' ${mcool}::/resolutions/1000 5000 10000 25000  \
                     --bigwig --min-dist-bad-bin 2

# Resolution 2,000 bp
cooltools insulation -o ${OUTDIR}/${species}_2000bpRes.boundaries_Li.tsv \
                     --threshold 'Li' ${mcool}::/resolutions/2000 10000 20000 50000  \
                     --bigwig --min-dist-bad-bin 2

# Resolution 4,000 bp
cooltools insulation -o ${OUTDIR}/${species}_4000bpRes.boundaries_Li.tsv \
                     --threshold 'Li' ${mcool}::/resolutions/4000 20000 40000 100000  \
                     --bigwig --min-dist-bad-bin 2