# early-metazoa-3D-chromatin
Code to reproduce analyses in Kim et al. 2024


### Data analyses
+ The **Micro-C processing pipeline**, which covers the full workflow from raw sequenced reads to normalized contact matrices, is available in the [scripts folder](./scripts/microc_processing/). This pipeline is adapted from the [4DN HiC processing pipeline](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline).
+ **Genome compartment analysis**, wihch can be found in the [script folder](./scripts/compartmentalization/) examines global compartmentalization patterns across phylogenetically diverse species. 
+ **Insulation analysis** includes the identification and annotation of insulation boundaries with epigenetic and genomic features. Find the relevant scripts [here](./scripts/insulation/)
+ The **Chromatin loop identification and annotation pipeline** for non-bilaterian animals is available [here](./scripts/chromatin_loops/).

### Data availability
+ Raw and processed files are available in GEO repository under accession number [GSE260572]().
+ All public datasets used in this study are listed [here](./data/Supplementary_Table_2_Public_datasets.xlsx).
+ _De novo_ sequenced and scaffolded to chromosome-level genomes of _Ephydatia muelleri_ and _Mnemiopsis leidyi_ are deposited in GenBank under [project number](). The genome of _Capsaspora owczarzaki_ can be found in DDBJ under [project number](). The same genomes are also availalbe in the [genome folder](./data/genome/).
+ Chromosome-level re-assemblies of are available in the [genome folder](./data/genome/).
+ Gene annotations and transposon annotations are deposited in the [data folder](./data/).
+ Identified insulation boundaries, genome compartments and chromatin loops are in the [data folder](./data/).

### Data exploration
Generated datasets could be explored in interactive genome browsers configured for each species:
[GBrowser](https://sebelab.crg.eu/3d-genomes-arc-jb2)