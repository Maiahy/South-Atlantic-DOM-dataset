# Data preprocessing

This folder contains the preprocessing workflows used for the South Atlantic DOM metabolomics dataset.

## Overview

Data preprocessing was performed in two independent steps:

1. **Full MS data preprocessing**, performed locally on a workstation.
2. **Data-dependent acquisition (DDA) data preprocessing**, performed on a high-performance computing cluster using SLURM.

The scripts provided in this folder are reproducible but require the preparation of the raw data and metadata folders as indicated in each script.

The raw LC-MS/MS data and associated metadata are publicly available through MetaboLights:

https://www.ebi.ac.uk/metabolights/MTBLS14508

## Output

At the end of this preprocessing stage, two independent feature tables are generated:

* a feature table corresponding to the Full MS dataset;
* a feature table corresponding to the DDA dataset.

These two datasets are subsequently linked and compared using the workflow available in the folder:

`2_Correspondance_DDA_fullMS`

## External function dependency

The preprocessing workflows require the file `customFunctions.R`, originally developed by Johannes Rainer as part of the xcms-gnps-tools project:

https://github.com/jorainer/xcms-gnps-tools

The version included in this repository was accessed on 2024-11-07 and is reproduced here unchanged in order to ensure long-term reproducibility and traceability of the computational workflow used in this study.

All credit for the development of `customFunctions.R` belongs to Johannes Rainer and contributors to the xcms-gnps-tools project. The file is included solely to facilitate reproducibility of the analyses presented in this repository.
