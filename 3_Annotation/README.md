# CAMERA annotation

This repository contains the code used to annotate isotopes and adducts in the Full MS dataset using the R package CAMERA.

The workflow starts from the preprocessed object generated during the data preprocessing step (`1_Data_preprocessing`) and produces a table containing CAMERA annotations.

Input file:

* `myEnvironment5_gap_filling_positive1.RData`

Output file:

* `result_CAMERA.csv`

Intermediate environments generated during the workflow are stored in `2_Environments`.

Raw data and metadata are available through MetaboLights:

https://www.ebi.ac.uk/metabolights/MTBLS14508

## Additional annotations

Annotations of the DDA dataset were performed independently using GNPS2 and Feature-Based Molecular Networking (FBMN).

The corresponding results are available at:

https://gnps2.org/status?task=1a38142f7fdb4cac919735abf489d06f

