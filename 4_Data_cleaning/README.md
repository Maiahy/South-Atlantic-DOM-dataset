# Data cleaning and normalization

This repository contains the code used to perform quality control, filtering, imputation, and batch correction of the Full MS dataset.

The workflow is divided into three scripts that must be run sequentially:

* `0A_First_steps_data_cleaning.R`
* `0B_App_SERRF.R`
* `0C_Last_part_data_cleaning.R`

The workflow starts from the full MS feature table generated during the XCMS preprocessing step and produces the final processed feature tables used for downstream statistical analyses.

Input files:

* Full MS feature table
* Sample metadata

Intermediate files generated during the workflow:

* `myEnvironment7_imputation2E5.RData`
* `positions_na`
* `normalized_by_SERRF.csv`

Output files:

* `0_final_feature_table_after_SERRF_RSD20.txt`
* `feature_table_SERRF_RSD20_with_NAs.txt`

Intermediate environments generated during the workflow are stored in the corresponding output folders.
