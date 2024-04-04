# early-metazoa-3D-chromatin
Code to reproduce analyses in Kim et al. 2024

### Analyses
+ Micro-C processing pipeline covering the steps from raw sequenced reads to normalized contact matrices is deposited [here](./scripts/microc_processing/). The pipeline is adapted from the [4DN HiC](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline) processing pipeline.
+ Genome [compartment analysis](./scripts/compartmentalization/). 
+ [Insulation analysis](./scripts/insulation/).
+ Identification and annotaiton of [chromatin loops](./scripts/chromatin_loops/).

### Data availability
+ Raw and processed files are available in GEO repository under accession number GSE260572.
+ All public datasets used in this study are listed [here](./data/Public_datasets_used_in_the_study.xlsx).
+ Genomes, gene annotations and transposon annotations are deposited in the [data folder](./data/).
+ Identified insulation boundaries, genome compartments and chromatin loops are in the [data folder](./data/).

### Data exploration
Generated datasets could be explored in interactive genome browsers configured for each species:
[Cowc_GBrowser](https://sebelab.crg.eu/cowc-hiroshi-jb2/) for _Capsaspora owczarzaki_<br>
[Emue_GBrowser](https://sebelab.crg.eu/emue-jb2/) for _Ephydatia muelleri_<br>
[Mlei_GBrowser](https://sebelab.crg.eu/mnemiopsis-jb2/) for _Mnemiopsis leidyi_<br>
[Tadh_GBrowser](https://sebelab.crg.eu/tadh-jb2/) for _Trichoplax adhaerens_
