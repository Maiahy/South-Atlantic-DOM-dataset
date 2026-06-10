# Correspondence between Full MS and DDA datasets

This folder contains the code used to link together the Full MS and DDA datasets using the R package **metabCombiner**.

The workflow follows the methodology described in the official metabCombiner vignette:

https://www.bioconductor.org/packages/release/bioc/vignettes/metabCombiner/inst/doc/metabCombiner_vignette.html

## Input files

The workflow requires two feature tables generated during the preprocessing step:

* `feature_table_MS1_positive1.txt`
* `DDAfeature_table_positive1.txt`

These files should be placed in the folder:

`1_Input_files`

The raw data and metadata used to generate these feature tables are available through MetaboLights:

https://www.ebi.ac.uk/metabolights/MTBLS14508

## Output files

Quality control figures are saved in:

`2_Quality_graphs`

Correspondence tables are saved in:

`3_Output_files`

## Notes

This workflow was used to establish correspondences between features detected in the Full MS and DDA datasets prior to downstream analyses.

