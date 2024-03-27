# **Compartmentalization analysis**

### To identify compartmentalization profiles, we followed the [cooltools notebook](https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html) workflow.

#### 1. Perform eigen value decomposition on balanced cooler matrices at resolutions equivalent to 5,000, 10,000, 20,000, 40,000 and 50,000 bins per species genome:
>
>[saddle_gc.ipynb](./saddle_gc.ipynb)
>

#### 2. Calculate compartment strength across different resolutions:
>
>[compartment_strength.R](./compartment_strength.R)
>

#### 3. Assign intermediate compartment to regions with weak compartmentalization:
>
>[intermediate_compartment.R](./intermediate_compartment.R)
>

#### 4. Calculate the distribution of genomic features, ChIP-seq chromatin signals and RNA-seq expression values across identified compartments:
>
>[omics_in_compartments.R](omics_in_compartments.R)
>
