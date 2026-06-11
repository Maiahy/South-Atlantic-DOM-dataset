# Technical replicate assessment

This repository contains the code used to evaluate the reproducibility of technical triplicate injections and to assess the influence of technical factors on sample distribution.

## Workflow

The workflow includes:

- Identification of samples showing high dispersion among technical triplicates using distances to sample centroids in Bray–Curtis PCoA space.
- Evaluation of the impact of removing highly dispersed samples.
- Comparison of PCoA structures before and after averaging technical triplicates.
- Assessment of the influence of evaporation status, in-lab operator, and sampling latitude on sample distribution.

## Scripts

### 1_Technical_replicates.R

Technical replicate assessment, outlier detection, and comparison of PCoA structures before and after filtering or averaging replicates.

### Operator A alone.R

Evaluation of technical and environmental covariates using only samples processed by operator A.

### Operator B alone.R

Evaluation of technical and environmental covariates using only samples processed by operator B.

## Input files

- `myEnvironment2_RSD_after_SERRF.RData`
- `feature_table_MS1_positive1.txt`

## Output files

- `PCoA_Distance_to_centroid.pdf`
- `Figure_SX_PCoA_processing.pdf`
- `Figure_SX_PCoA_covariantsAf.pdf`
- `Figure_SX_PCoA_covariantsBf.pdf`

## Data availability

Raw data and metadata are available through MetaboLights:

https://www.ebi.ac.uk/metabolights/MTBLS14508
