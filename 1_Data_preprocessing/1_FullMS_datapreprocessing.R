
#Full MS dataset 

#this code is based on the work of different groups, for more information please consult:
#https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html
#https://github.com/DorresteinLaboratory/XCMS3_FeatureBasedMN/blob/master/XCMS3_Preprocessing.Rmd
#https://pietrofranceschi.github.io/Metabolomics_lectures/
#https://rformassspectrometry.github.io/Metabonaut/articles/a-end-to-end-untargeted-metabolomics.html

#####1) Packages loading ####

#install.packages("BiocManager")
#BiocManager::install("MsExperiment")
#dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
#.libPaths(Sys.getenv("R_LIBS_USER"))
#install.packages("remotes")
#remotes::install_github("RforMassSpectrometry/MsIO")
#BiocManager::install("alabaster.se")
#BiocManager::install("RforMassSpectrometry/MsBackendMetaboLights")
#BiocManager::install("MsBackendMgf")
#BiocManager::install("xcms")
#BiocManager::install("AnnotationHub")
#BiocManager::install("CompoundDb")
#BiocManager::install("MetaboAnnotation")
#remotes::install_github("https://github.com/girke-lab/ChemmineR") # You need the last version of Rtool on your path + do not download the libraries or it won't work fo some reasons...
#install.packages("janitor")


## Data Import and handling
library(knitr)
library(readxl)
library(MsExperiment)
library(MsIO)
library(alabaster.se)
library(MsBackendMetaboLights)
library(SummarizedExperiment)
library(MsBackendMgf)
library(janitor)
library(dplyr)

## Preprocessing of LC-MS data
library(xcms)
library(Spectra)
library(MetaboCoreUtils)

## Statistical analysis
library(limma) # Differential abundance
library(matrixStats) # Summaries over matrices

## Visualisation
library(pander)
library(RColorBrewer)
library(pheatmap)
library(vioplot)
library(ggfortify) # Plot PCA
library(gridExtra) # To arrange multiple ggplots into single plots

## Annotation
library(ChemmineR)
library(CompoundDb) # Access small compound annotation data.
library(MetaboAnnotation) # Functionality for metabolite annotation.
library(AnnotationHub)# Annotation resources

####2) FULL MS  ####
#####2.a) Paths ####

#Define main path
library(here)
base_path <- here()

#Define sub-paths
raw_data_path <- file.path(base_path, "1_Raw_data")
metadata_path <- file.path(base_path, "Metadata_final_with_operators.txt")
out_path <- file.path(base_path, "2_Output_files")
env_path <- file.path(out_path, "1_Environments")
plot_path <- file.path(out_path, "2_Quality_plots")

#Create directories (you need to have 1_Raw_data and Metadata_final_with_operators.txt created)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(env_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

#Select a polarity:
polarity_mode <- 1  # 1 for positive mode, 0 for negative mode.

#####2.b) Data import full MS (mse1) ####
#Read metadata 
metadata1 <- read.table(metadata_path, sep = "\t") %>%
  row_to_names(row_number = 1) 
colnames(metadata1)
unique(metadata1$batch_identifier)

#List files in the raw data directory
mydata1 <- list.files(path = raw_data_path, full.names = TRUE) 

#Remove the files from DDA dataset
exclusion_patterns <- c("_MALLMB02B1POS_", "_MALLMB02B2POS_", "_MALLMB02B3POS_", 
                        "_MALLMB02B4POS_","_B1_", "_B2_", "_B3_", "_B4_", "_B5_")
exclusion_pattern <- paste(exclusion_patterns, collapse = "|")
filtered_files <- mydata1[!grepl(exclusion_pattern, mydata1, ignore.case = TRUE)]

#Create a MsExperiment object
pd <- data.frame(file_name = basename(filtered_files))
pd <- left_join(pd, metadata1, by = "file_name")
mse1 <- readMsExperiment(filtered_files, sampleData = pd)

#Filter raw data to keep only the selected polarity
mydata1 <- filterSpectra(mse1, filterPolarity, polarity = polarity_mode)
spectra_data <- spectra(mydata1)
ms_polarity_table <- table(msLevel(spectra_data), polarity(spectra_data))
ms_polarity_df <- as.data.frame(ms_polarity_table)
write.csv(ms_polarity_df, file = file.path(out_path, paste0("msLevel_vs_polarity_postive1", ".csv")), 
          row.names = FALSE)

save.image(
  file = file.path(
    env_path,
    "myEnvironment1_Data_Loading1.RData"
  )
)

#####2.c) Data organisation ####

#Open and analyze your object:
mydata1

#Object of class MsExperiment 
#Spectra: MS1 (753803) 
#Experiment data: 1393 sample(s)
#Sample data links:
#  - spectra: 1393 sample(s) to 753803 element(s).

sampleData(mydata1) #the metadata are stored in "sampleData"
spectra(mydata1) #the metadata are stored in "spectra"

#Extract and plot BPC for full data
bpc <- chromatogram(mse1, aggregationFun = "max")
plot(bpc)

register(SerialParam())
set.seed(070499)
random_indices <- sample(seq_len(length(mydata1)), 50)
subset_mydata1 = mydata1[random_indices]
get_mz_window <- function(mz_center, ppm = 100) {
  return((ppm * mz_center) / 1e6)
}

#Liste of the priority compounds
priority_compounds <- list(
  list(name = "Choline", mz = 104.11, rt = c(0, 200)),
  list(name = "Betaine", mz = 118.09, rt = c(60.21 - 120, 60.21 + 120)),
  list(name = "Homarine", mz = 138.06, rt = c(60.21 - 120, 60.21 + 120)),
  list(name = "Salicylic_acid", mz = 121.03, rt = c(346.38 - 120, 346.38 + 120)),
  list(name = "Nicotinic_acid", mz = 124.04, rt = c(66.45 - 120, 66.45 + 120)),
  list(name = "Arginine", mz = 175.12, rt = c(53.35 - 120, 53.35 + 120)),
  list(name = "Val-Leu", mz = 231.17, rt = c(183.02 - 120, 183.02 + 120)),
  list(name = "Thr_Leu", mz = 233.15, rt = c(170.58 - 120, 170.58 + 120)),
  list(name = "Stachydrine", mz = 144.10, rt = c(60.84 - 120, 60.84 + 120)),
  list(name = "Tyrosine", mz = 182.08, rt = c(144.09 - 120, 144.09 + 120)),
  list(name = "Caffeine", mz = 195.09, rt = c(199.16 - 120, 199.16 + 120)),
  list(name = "Vanillin", mz = 153.06, rt = c(243.39 - 120, 243.39 + 120)),
  list(name = "homocysteine", mz = 385.13, rt = c(72.68 - 120, 72.68 + 120)),
  list(name = "Guanine", mz = 152.06, rt = c(110.74 - 120, 110.74 + 120)),
  list(name = "Loliolide", mz = 197.12, rt = c(259.28 - 120, 259.28 + 120)),
  list(name = "Myrcene", mz = 159.12, rt = c(388.85 - 120, 388.85 + 120)),
  list(name = "Astaxanthin", mz = 597.39, rt = c(606.32 - 120, 606.32 + 120)),
  list(name = "Lumichrome", mz = 243.09, rt = c(256.55 - 120, 256.55 + 120)),
  list(name = "Diadinoxanthin", mz = 581.40, rt = c(567.38 - 120, 567.38 + 120)),
  list(name = "Crassostreaxanthin", mz = 599.41, rt = c(563.04 - 120, 563.04 + 120)),
  list(name = "Shinorine", mz = 333.13, rt = c(62.44 - 120, 62.44 + 120)),
  list(name = "Pinolenic_acid", mz = 161.10, rt = c(277.36 - 120, 277.36 + 120)),
  list(name = "Trimethyllysine", mz = 189.16, rt = c(55.49 - 120, 55.49 + 120)),
  list(name = "Usujirene", mz = 285.14, rt = c(148.69 - 120, 148.69 + 120)),
  list(name = "Anthranilate", mz = 120.04, rt = c(81.32 - 120, 81.32 + 120)),
  list(name = "Succinoadenosine", mz = 384.12, rt = c(168.22 - 120, 168.22 + 120)),
  list(name = "Ala-Phe", mz = 237.12, rt = c(183.02 - 120, 183.02 + 120)),
  list(name = "Diphenylphosphate", mz = 251.05, rt = c(521.71 - 120, 521.71 + 120)),
  list(name = "Octinoxate", mz = 291.19, rt = c(438.85 - 120, 438.85 + 120)),
  list(name = "Val-Trp", mz = 304.17, rt = c(209.34 - 120, 209.34 + 120)),
  list(name = "Inosine", mz = 269.09, rt = c(92.87 - 120, 92.87 + 120)),
  list(name = "Methyl_coumarin", mz = 161.06, rt = c(317.59 - 120, 317.59 + 120)),
  list(name = "Smenamide", mz = 501.25, rt = c(432.08 - 120, 432.08 + 120)),
  list(name = "Ser-Ile", mz = 219.13, rt = c(153.19 - 120, 153.19 + 120)),
  list(name = "Phe-Glu", mz = 295.13, rt = c(161.04 - 120, 161.04 + 120)),
  list(name = "Petasol", mz = 235.17, rt = c(453.49 - 120, 453.49 + 120)),
  list(name = "Butyrylcarnitine", mz = 232.15, rt = c(168.22 - 120, 168.22 + 120)),
  list(name = "palmitoyl_carnitine", mz = 400.34, rt = c(476.56 - 120, 476.56 + 120)),
  list(name = "Benzothiazolesulfonic_acid", mz = 215.98, rt = c(211.42 - 120, 211.42 + 120)),
  list(name = "LAUROYLCARNITINE", mz = 344.28, rt = c(369.72 - 120, 369.72 + 120)),
  list(name = "Phosphocholine", mz = 184.07, rt = c(397.18 - 120, 397.18 + 120))
)


#EIC plot
for (compound in priority_compounds) {
  eic <- chromatogram(subset_mydata1,
                      mz = compound$mz + c(-get_mz_window(compound$mz), get_mz_window(compound$mz)),
                      rt = compound$rt)
  plot(eic, main = paste("EIC for", compound$name, "(±100 ppm)"),
       sub = paste("m/z =", compound$mz, "| RT =", compound$rt[1], "-", compound$rt[2]))
}



#####2.d) Peak picking ####

#Set up parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(MulticoreParam(2))
} else {
  register(SnowParam(2))
} #In case you can use parallel processing :)

register(SerialParam()) #In case you can't for some reasons

#Set the parameters
#noise = 2E5 estimated with mzmine3 (v3.2.8) 
#peakwidth optimised according the same method
cwp <- CentWaveParam(snthresh = 5, noise = 2E5, peakwidth = c(2, 20),
                     ppm = 15, prefilter = c(3, 3000), integrate = 2)

#Perform the chromatographic peak detection using *centWave*.
mse1 <- findChromPeaks(mydata1, param = cwp)


save.image(
  file = file.path(
    env_path,
    "myEnvironment2_Peak_picking_positive1.RData"
  )
)

chromPeaks(mse1) |>
  head()

register(SerialParam())
#Check for some compounds of interest:

#Choline?
eic_choline <- chromatogram(mse1, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70))
#Adenosine?
eic_Adenosine <- chromatogram(mse1, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60))
#Phosphocholine?
eic_Phosphocholine <- chromatogram(mse1, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60))
#Nicotinic_acid?
eic_Nicotinic_acid <- chromatogram(mse1, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60))
#Arginine
eic_Arginine <- chromatogram(mse1, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60))
#Val-Leu
eic_Val_Leu <- chromatogram(mse1, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60))
#Homarine
eic_Homarine <- chromatogram(mse1, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60))
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(mse1, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60))
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(mse1, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60))
#Diphenylamine
eic_Diphenylamine <- chromatogram(mse1, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60))
#Thr_Leu
eic_Thr_Leu <- chromatogram(mse1, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60))

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(mse1, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60))
#Betain_lipid
eic_Betain_lipid <- chromatogram(mse1, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60))

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(mse1, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60))

#Tyrosine
eic_Tyrosine <- chromatogram(mse1, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60))

#Caffeine
eic_Caffeine <- chromatogram(mse1, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60))
#Vanillin
eic_Vanillin <- chromatogram(mse1, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60))
#Stachydrine
eic_Stachydrine <- chromatogram(mse1, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60))
#homocysteine
eic_homocysteine <- chromatogram(mse1, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-60, 60))
#Guanine
eic_Guanine <- chromatogram(mse1, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60))
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(mse1, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60))
#Riboflavin
eic_Riboflavin <- chromatogram(mse1, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60))
#Cholesterol
eic_Cholesterol <- chromatogram(mse1, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60))
#Ile_Lys
eic_Ile_Lys <- chromatogram(mse1, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60))
#Loliolide
eic_Loliolide <- chromatogram(mse1, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60))
#myrcene
eic_Myrcene <- chromatogram(mse1, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60))
#Astaxanthin
eic_Astaxanthin <- chromatogram(mse1, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60))
#LUMICHROME
eic_Lumichrome <- chromatogram(mse1, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60))
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(mse1, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60))
#BETAINE
eic_Betain <- chromatogram(mse1, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60))
#BiakamideC,D
eic_Biakamide <- chromatogram(mse1, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60))
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(mse1, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60))
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(mse1, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60))
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(mse1, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60))
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(mse1, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60))
#Shinorine
eic_Shinorine  <- chromatogram(mse1, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60))
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(mse1, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60))
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(mse1, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60))
#Usujirene
eic_Usujirene <- chromatogram(mse1, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60))
#ANTHRANILATE
eic_Anthranilate <- chromatogram(mse1, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60))
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(mse1, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#Ala-Phe
eic_Ala_Phe <- chromatogram(mse1, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60))
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(mse1, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60))
#Octinoxate
eic_Octinoxate <- chromatogram(mse1, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60))
#Val-Trp
eic_Val_Trp <- chromatogram(mse1, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60))
#Inosine
eic_Inosine <- chromatogram(mse1, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60))
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(mse1, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60))
#Smenamide A
eic_Smenamide <- chromatogram(mse1, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60))
#Ser-Ile
eic_Ser_Ile <- chromatogram(mse1, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60))
#Phe-Glu
eic_Phe_Glu <- chromatogram(mse1, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60))
#Petasol
eic_Petasol <- chromatogram(mse1, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60))
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(mse1, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(mse1, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60))
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(mse1, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60))
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(mse1, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60))

save.image(
  file = file.path(
    env_path,
    "myEnvironment2_Peak_picking_EIC_positive1.RData"
  )
)

library(grDevices)
pdf(
  file.path(
    plot_path,
    "Peak_picking_quality_positive1.pdf"
  ),
  width = 14,
  height = 10
)
par(mfrow = c(4, 5))
peak_bg_color <- adjustcolor("cyan4", alpha.f = 0.3)
# Plot each compound with blue peak areas
plot(eic_choline, main = expression(paste(italic(m/z), " 104.11")), peakBg = peak_bg_color)
plot(eic_Adenosine, main = expression(paste(italic(m/z), " 268.10")), peakBg = peak_bg_color)
plot(eic_Phosphocholine, main = expression(paste(italic(m/z), " 184.07")), peakBg = peak_bg_color)
plot(eic_Nicotinic_acid, main = expression(paste(italic(m/z), " 124.04")), peakBg = peak_bg_color)
plot(eic_Arginine, main = expression(paste(italic(m/z), " 175.12")), peakBg = peak_bg_color)
plot(eic_Val_Leu, main = expression(paste(italic(m/z), " 231.17")), peakBg = peak_bg_color)
plot(eic_Homarine, main = expression(paste(italic(m/z), " 138.06")), peakBg = peak_bg_color)
plot(eic_Salicylic_acid, main = expression(paste(italic(m/z), " 121.03")), peakBg = peak_bg_color)
plot(eic_Lauryl_diethanolamide, main = expression(paste(italic(m/z), " 288.25")), peakBg = peak_bg_color)
plot(eic_Diphenylamine, main = expression(paste(italic(m/z), " 170.10")), peakBg = peak_bg_color)
plot(eic_Thr_Leu, main = expression(paste(italic(m/z), " 233.15")), peakBg = peak_bg_color)
plot(eic_CocamidoprpylBetaine, main = expression(paste(italic(m/z), " 343.30")), peakBg = peak_bg_color)
plot(eic_Betain_lipid, main = expression(paste(italic(m/z), " 490.37")), peakBg = peak_bg_color)
plot(eic_p_coumaric_acid, main = expression(paste(italic(m/z), " 147.04")), peakBg = peak_bg_color)
plot(eic_Tyrosine, main = expression(paste(italic(m/z), " 182.08")), peakBg = peak_bg_color)
plot(eic_Caffeine, main = expression(paste(italic(m/z), " 195.09")), peakBg = peak_bg_color)
plot(eic_Vanillin, main = expression(paste(italic(m/z), " 153.06")), peakBg = peak_bg_color)
plot(eic_Stachydrine, main = expression(paste(italic(m/z), " 144.10")), peakBg = peak_bg_color)
plot(eic_homocysteine, main = expression(paste(italic(m/z), " 385.13")), peakBg = peak_bg_color)
plot(eic_Guanine, main = expression(paste(italic(m/z), " 152.06")), peakBg = peak_bg_color)
plot(eic_Dinoxanthin, main = expression(paste(italic(m/z), " 641.42")), peakBg = peak_bg_color)
plot(eic_Riboflavin, main = expression(paste(italic(m/z), " 377.15")), peakBg = peak_bg_color)
plot(eic_Cholesterol, main = expression(paste(italic(m/z), " 369.35")), peakBg = peak_bg_color)
plot(eic_Ile_Lys, main = expression(paste(italic(m/z), " 260.20")), peakBg = peak_bg_color)
plot(eic_Loliolide, main = expression(paste(italic(m/z), " 197.12")), peakBg = peak_bg_color)
plot(eic_Myrcene, main = expression(paste(italic(m/z), " 159.12")), peakBg = peak_bg_color)
plot(eic_Astaxanthin, main = expression(paste(italic(m/z), " 597.39")), peakBg = peak_bg_color)
plot(eic_Lumichrome, main = expression(paste(italic(m/z), " 243.09")), peakBg = peak_bg_color)
plot(eic_Diadinoxanthin, main = expression(paste(italic(m/z), " 581.40")), peakBg = peak_bg_color)
plot(eic_Betain, main = expression(paste(italic(m/z), " 118.09")), peakBg = peak_bg_color)
plot(eic_Biakamide, main = expression(paste(italic(m/z), " 524.27")), peakBg = peak_bg_color)
plot(eic_Fucoxanthin, main = expression(paste(italic(m/z), " 659.43")), peakBg = peak_bg_color)
plot(eic_Arachidonic_acid, main = expression(paste(italic(m/z), " 305.25")), peakBg = peak_bg_color)
plot(eic_Phaeophorbide , main = expression(paste(italic(m/z), " 593.27")), peakBg = peak_bg_color)
plot(eic_Crassostreaxanthin , main = expression(paste(italic(m/z), " 599.41")), peakBg = peak_bg_color)
plot(eic_Shinorine , main = expression(paste(italic(m/z), " 333.13")), peakBg = peak_bg_color)
plot(eic_Pinolenic_acid , main = expression(paste(italic(m/z), " 161.10")), peakBg = peak_bg_color)
plot(eic_Trimethyllysine , main = expression(paste(italic(m/z), " 189.16")), peakBg = peak_bg_color)
plot(eic_Usujirene , main = expression(paste(italic(m/z), " 285.14")), peakBg = peak_bg_color)
plot(eic_Anthranilate , main = expression(paste(italic(m/z), " 120.04")), peakBg = peak_bg_color)
plot(eic_Succinoadenosine , main = expression(paste(italic(m/z), " 384.12")), peakBg = peak_bg_color)
plot(eic_Ala_Phe , main = expression(paste(italic(m/z), " 237.12")), peakBg = peak_bg_color)
plot(eic_Diphenylphosphate , main = expression(paste(italic(m/z), " 251.05")), peakBg = peak_bg_color)
plot(eic_Octinoxate , main = expression(paste(italic(m/z), " 291.19")), peakBg = peak_bg_color)
plot(eic_Val_Trp , main = expression(paste(italic(m/z), " 304.17")), peakBg = peak_bg_color)
plot(eic_Inosine , main = expression(paste(italic(m/z), " 269.09")), peakBg = peak_bg_color)
plot(eic_Methyl_coumarin , main = expression(paste(italic(m/z), " 161.06")), peakBg = peak_bg_color)
plot(eic_Smenamide , main = expression(paste(italic(m/z), " 501.25")), peakBg = peak_bg_color)
plot(eic_Ser_Ile , main = expression(paste(italic(m/z), " 219.13")), peakBg = peak_bg_color)
plot(eic_Phe_Glu , main = expression(paste(italic(m/z), " 295.13")), peakBg = peak_bg_color)
plot(eic_Petasol , main = expression(paste(italic(m/z), " 235.17")), peakBg = peak_bg_color)
plot(eic_Butyrylcarnitine , main = expression(paste(italic(m/z), " 232.15")), peakBg = peak_bg_color)
plot(eic_palmitoyl_carnitine , main = expression(paste(italic(m/z), " 400.34")), peakBg = peak_bg_color)
plot(eic_Benzothiazolesulfonic_acid , main = expression(paste(italic(m/z), " 215.98")), peakBg = peak_bg_color)
plot(eic_Lauroylcarnitine , main = expression(paste(italic(m/z), " 344.28")), peakBg = peak_bg_color)
dev.off()

#####2.e) Peak refinement ####

#Parallel processing using 2 cores
if (.Platform$OS.type == "unix") {
  register(MulticoreParam(2))
} else {
  register(SnowParam(2))
} #(In case you can use parallel processing!)

#register(SerialParam()) #In case you can't for some reasons

param <- MergeNeighboringPeaksParam(expandRt = 2.5, expandMz = 0.0015,
                                    minProp = 0.75)
mse1 <- refineChromPeaks(mse1, param = param, chunkSize = 5)

save.image(
  file = file.path(
    env_path,
    "myEnvironment2_Peak_refinement_positive1.RData"
  )
)

register(SerialParam()) 
#Choline?
eic_choline <- chromatogram(mse1, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70))
#Adenosine?
eic_Adenosine <- chromatogram(mse1, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60))
#Phosphocholine?
eic_Phosphocholine <- chromatogram(mse1, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60))
#Nicotinic_acid?
eic_Nicotinic_acid <- chromatogram(mse1, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60))
#Arginine
eic_Arginine <- chromatogram(mse1, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60))
#Val-Leu
eic_Val_Leu <- chromatogram(mse1, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60))
#Homarine
eic_Homarine <- chromatogram(mse1, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60))
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(mse1, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60))
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(mse1, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60))
#Diphenylamine
eic_Diphenylamine <- chromatogram(mse1, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60))
#Thr_Leu
eic_Thr_Leu <- chromatogram(mse1, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60))

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(mse1, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60))
#Betain_lipid
eic_Betain_lipid <- chromatogram(mse1, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60))

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(mse1, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60))

#Tyrosine
eic_Tyrosine <- chromatogram(mse1, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60))

#Caffeine
eic_Caffeine <- chromatogram(mse1, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60))
#Vanillin
eic_Vanillin <- chromatogram(mse1, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60))
#Stachydrine
eic_Stachydrine <- chromatogram(mse1, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60))
#homocysteine
eic_homocysteine <- chromatogram(mse1, mz = 	385.13 + c(-0.01, 0.01),
                               rt = 72.68 + c(-60, 60))
#Guanine
eic_Guanine <- chromatogram(mse1, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60))
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(mse1, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60))
#Riboflavin
eic_Riboflavin <- chromatogram(mse1, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60))
#Cholesterol
eic_Cholesterol <- chromatogram(mse1, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60))
#Ile_Lys
eic_Ile_Lys <- chromatogram(mse1, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60))
#Loliolide
eic_Loliolide <- chromatogram(mse1, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60))
#myrcene
eic_Myrcene <- chromatogram(mse1, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60))
#Astaxanthin
eic_Astaxanthin <- chromatogram(mse1, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60))
#LUMICHROME
eic_Lumichrome <- chromatogram(mse1, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60))
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(mse1, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60))
#BETAINE
eic_Betain <- chromatogram(mse1, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60))
#BiakamideC,D
eic_Biakamide <- chromatogram(mse1, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60))
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(mse1, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60))
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(mse1, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60))
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(mse1, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60))
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(mse1, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60))
#Shinorine
eic_Shinorine  <- chromatogram(mse1, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60))
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(mse1, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60))
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(mse1, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60))
#Usujirene
eic_Usujirene <- chromatogram(mse1, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60))
#ANTHRANILATE
eic_Anthranilate <- chromatogram(mse1, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60))
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(mse1, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#Ala-Phe
eic_Ala_Phe <- chromatogram(mse1, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60))
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(mse1, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60))
#Octinoxate
eic_Octinoxate <- chromatogram(mse1, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60))
#Val-Trp
eic_Val_Trp <- chromatogram(mse1, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60))
#Inosine
eic_Inosine <- chromatogram(mse1, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60))
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(mse1, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60))
#Smenamide A
eic_Smenamide <- chromatogram(mse1, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60))
#Ser-Ile
eic_Ser_Ile <- chromatogram(mse1, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60))
#Phe-Glu
eic_Phe_Glu <- chromatogram(mse1, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60))
#Petasol
eic_Petasol <- chromatogram(mse1, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60))
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(mse1, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(mse1, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60))
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(mse1, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60))
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(mse1, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60))

save.image(
  file = file.path(
    env_path,
    "myEnvironment2_Peak_refinement_EIC_positive1.RData"
  )
)

library(grDevices)
pdf(
  file.path(
    plot_path,
    "Peak_refinement_quality_positive1.pdf"
  ),
  width = 14,
  height = 10
)
par(mfrow = c(4, 5))
peak_bg_color <- adjustcolor("cyan4", alpha.f = 0.3)

# Plot each compound with blue peak areas
plot(eic_choline, main = expression(paste(italic(m/z), " 104.11")), peakBg = peak_bg_color)
plot(eic_Adenosine, main = expression(paste(italic(m/z), " 268.10")), peakBg = peak_bg_color)
plot(eic_Phosphocholine, main = expression(paste(italic(m/z), " 184.07")), peakBg = peak_bg_color)
plot(eic_Nicotinic_acid, main = expression(paste(italic(m/z), " 124.04")), peakBg = peak_bg_color)
plot(eic_Arginine, main = expression(paste(italic(m/z), " 175.12")), peakBg = peak_bg_color)
plot(eic_Val_Leu, main = expression(paste(italic(m/z), " 231.17")), peakBg = peak_bg_color)
plot(eic_Homarine, main = expression(paste(italic(m/z), " 138.06")), peakBg = peak_bg_color)
plot(eic_Salicylic_acid, main = expression(paste(italic(m/z), " 121.03")), peakBg = peak_bg_color)
plot(eic_Lauryl_diethanolamide, main = expression(paste(italic(m/z), " 288.25")), peakBg = peak_bg_color)
plot(eic_Diphenylamine, main = expression(paste(italic(m/z), " 170.10")), peakBg = peak_bg_color)
plot(eic_Thr_Leu, main = expression(paste(italic(m/z), " 233.15")), peakBg = peak_bg_color)
plot(eic_CocamidoprpylBetaine, main = expression(paste(italic(m/z), " 343.30")), peakBg = peak_bg_color)
plot(eic_Betain_lipid, main = expression(paste(italic(m/z), " 490.37")), peakBg = peak_bg_color)
plot(eic_p_coumaric_acid, main = expression(paste(italic(m/z), " 147.04")), peakBg = peak_bg_color)
plot(eic_Tyrosine, main = expression(paste(italic(m/z), " 182.08")), peakBg = peak_bg_color)
plot(eic_Caffeine, main = expression(paste(italic(m/z), " 195.09")), peakBg = peak_bg_color)
plot(eic_Vanillin, main = expression(paste(italic(m/z), " 153.06")), peakBg = peak_bg_color)
plot(eic_Stachydrine, main = expression(paste(italic(m/z), " 144.10")), peakBg = peak_bg_color)
plot(eic_homocysteine, main = expression(paste(italic(m/z), " 385.13")), peakBg = peak_bg_color)
plot(eic_Guanine, main = expression(paste(italic(m/z), " 152.06")), peakBg = peak_bg_color)
plot(eic_Dinoxanthin, main = expression(paste(italic(m/z), " 641.42")), peakBg = peak_bg_color)
plot(eic_Riboflavin, main = expression(paste(italic(m/z), " 377.15")), peakBg = peak_bg_color)
plot(eic_Cholesterol, main = expression(paste(italic(m/z), " 369.35")), peakBg = peak_bg_color)
plot(eic_Ile_Lys, main = expression(paste(italic(m/z), " 260.20")), peakBg = peak_bg_color)
plot(eic_Loliolide, main = expression(paste(italic(m/z), " 197.12")), peakBg = peak_bg_color)
plot(eic_Myrcene, main = expression(paste(italic(m/z), " 159.12")), peakBg = peak_bg_color)
plot(eic_Astaxanthin, main = expression(paste(italic(m/z), " 597.39")), peakBg = peak_bg_color)
plot(eic_Lumichrome, main = expression(paste(italic(m/z), " 243.09")), peakBg = peak_bg_color)
plot(eic_Diadinoxanthin, main = expression(paste(italic(m/z), " 581.40")), peakBg = peak_bg_color)
plot(eic_Betain, main = expression(paste(italic(m/z), " 118.09")), peakBg = peak_bg_color)
plot(eic_Biakamide, main = expression(paste(italic(m/z), " 524.27")), peakBg = peak_bg_color)
plot(eic_Fucoxanthin, main = expression(paste(italic(m/z), " 659.43")), peakBg = peak_bg_color)
plot(eic_Arachidonic_acid, main = expression(paste(italic(m/z), " 305.25")), peakBg = peak_bg_color)
plot(eic_Phaeophorbide , main = expression(paste(italic(m/z), " 593.27")), peakBg = peak_bg_color)
plot(eic_Crassostreaxanthin , main = expression(paste(italic(m/z), " 599.41")), peakBg = peak_bg_color)
plot(eic_Shinorine , main = expression(paste(italic(m/z), " 333.13")), peakBg = peak_bg_color)
plot(eic_Pinolenic_acid , main = expression(paste(italic(m/z), " 161.10")), peakBg = peak_bg_color)
plot(eic_Trimethyllysine , main = expression(paste(italic(m/z), " 189.16")), peakBg = peak_bg_color)
plot(eic_Usujirene , main = expression(paste(italic(m/z), " 285.14")), peakBg = peak_bg_color)
plot(eic_Anthranilate , main = expression(paste(italic(m/z), " 120.04")), peakBg = peak_bg_color)
plot(eic_Succinoadenosine , main = expression(paste(italic(m/z), " 384.12")), peakBg = peak_bg_color)
plot(eic_Ala_Phe , main = expression(paste(italic(m/z), " 237.12")), peakBg = peak_bg_color)
plot(eic_Diphenylphosphate , main = expression(paste(italic(m/z), " 251.05")), peakBg = peak_bg_color)
plot(eic_Octinoxate , main = expression(paste(italic(m/z), " 291.19")), peakBg = peak_bg_color)
plot(eic_Val_Trp , main = expression(paste(italic(m/z), " 304.17")), peakBg = peak_bg_color)
plot(eic_Inosine , main = expression(paste(italic(m/z), " 269.09")), peakBg = peak_bg_color)
plot(eic_Methyl_coumarin , main = expression(paste(italic(m/z), " 161.06")), peakBg = peak_bg_color)
plot(eic_Smenamide , main = expression(paste(italic(m/z), " 501.25")), peakBg = peak_bg_color)
plot(eic_Ser_Ile , main = expression(paste(italic(m/z), " 219.13")), peakBg = peak_bg_color)
plot(eic_Phe_Glu , main = expression(paste(italic(m/z), " 295.13")), peakBg = peak_bg_color)
plot(eic_Petasol , main = expression(paste(italic(m/z), " 235.17")), peakBg = peak_bg_color)
plot(eic_Butyrylcarnitine , main = expression(paste(italic(m/z), " 232.15")), peakBg = peak_bg_color)
plot(eic_palmitoyl_carnitine , main = expression(paste(italic(m/z), " 400.34")), peakBg = peak_bg_color)
plot(eic_Benzothiazolesulfonic_acid , main = expression(paste(italic(m/z), " 215.98")), peakBg = peak_bg_color)
plot(eic_Lauroylcarnitine , main = expression(paste(italic(m/z), " 344.28")), peakBg = peak_bg_color)
dev.off()

#####2.f) Retention time alignment ####

register(SerialParam())
mse0 <- mse1  # Backup before alignment for future QC figures

#Metadata inspection
metadata_check <- data.frame(
  file_name = sampleData(mse1)$file_name,
  ESSENTIAL_sample_id = sampleData(mse1)$ESSENTIAL_sample_id,
  ESSENTIAL_sample_type = sampleData(mse1)$ESSENTIAL_sample_type,
  stringsAsFactors = FALSE)
write.csv(metadata_check,
          file = file.path(out_path, "metadata_check_sampleType_positive1.csv"),
          row.names = FALSE, quote = TRUE)

#Sort dataset by injection order
sample_data <- sampleData(mse1)
sample_data$Injection_order <- as.numeric(sample_data$Injection_order)
sample_data <- sample_data[order(sample_data$Injection_order), ]
sampleData(mse1) <- sample_data  # Reassign sorted metadata

# Check injection order
print(sampleData(mse1)$Injection_order)

#Qulity insurance: injection order
injection_order_df <- data.frame(
  file_name = sample_data$file_name,
  injection_order = sample_data$Injection_order,
  ESSENTIAL_sample_type = sample_data$ESSENTIAL_sample_type)
write.csv(injection_order_df, file = file.path(out_path, paste0("Injection_order_check_positive1", ".csv")), 
          row.names = FALSE)

#Grouping
mse1 <- groupChromPeaks(mse1, param = PeakDensityParam(
    sampleGroups = sampleData(mse1)$ESSENTIAL_sample_id,
    minFraction = 2/3,           
    binSize = 0.02, #Rule: binSize >= max(mz * ppm / 1e6)
    ppm = 10,
    bw = 3))

#Adjust RT
mse1 <- adjustRtime(mse1, param = PeakGroupsParam(
                               minFraction = 0.8,
                               subset = which(sampleData(mse1)$ESSENTIAL_sample_id == "InterbatchQC"), 
                               subsetAdjust = "average",
                               span = 0.6))

# Save environment
save.image(
  file = file.path(
    env_path,
    "myEnvironment3_RT_alignement_positive1.RData"
  )
)
#Quality control figures:

#1)RT adjustment:
pd <- sampleData(mse1)
sample_colors <- c("QC" = "#213146","Sample" = "#579ad1")
color_mapping <- sample_colors[sampleData(mse1)$ESSENTIAL_sample_type]
pdf(
  file.path(
    plot_path,
    "Retention_time_alignment_positive_simplified.pdf"
  ),
  width = 14,
  height = 10
)
plotAdjustedRtime(mse1, col = color_mapping, main = "Adjusted retention times across samples")
legend("topright", legend = names(sample_colors), col = sample_colors, pch = 15, cex = 0.8, title = "Sample Type", bg = "white")
dev.off()


library(grDevices)

pd <- sampleData(mse1)

sample_type <- pd$ESSENTIAL_sample_type

# Create per-sample colors
color_mapping <- ifelse(
  sample_type == "QC",
  adjustcolor("#213146", alpha.f = 0.9),   # strong QC
  ifelse(
    sample_type == "Blank",
    adjustcolor("#df6765", alpha.f = 0.50),
    adjustcolor("#579ad1", alpha.f = 0.20)
  )
)

sample_colors <- c("QC" = "#213146","Sample" = "#579ad1", "Blank" = "#df6765")

#2)BPCs before and after alignment all samples:
bpc_raw <- chromatogram(mse0, aggregationFun = "max")
bpc_aligned <- chromatogram(mse1, aggregationFun = "max", include = "none")
pdf(
  file.path(
    plot_path,
    "Retention_time_alignment2_positivemap2.pdf"
  ),
  width = 14,
  height = 10
)
par(mfrow = c(2, 1))
plot(bpc_raw, peakType = "none", col = color_mapping, main = "Base Peak Chromatogram before Retention Time Alignment")
grid()
legend("topright", legend = names(sample_colors), col = sample_colors, pch = 15, cex = 0.8, title = "Sample Type", bg = "white")
plot(bpc_aligned, col = color_mapping, main = "Base Peak Chromatogram after Retention Time Alignment")
grid()
dev.off()

#2)Without the QCs:
#Only non-QC samples
keep <- sample_type != "QC"

#Subset chromatograms
bpc_raw <- chromatogram(mse0[keep], aggregationFun = "max")
bpc_aligned <- chromatogram(mse1[keep], aggregationFun = "max", include = "none")

#Colors without QC
color_mapping <- ifelse(
  sample_type[keep] == "Blank",
  adjustcolor("#df6765", alpha.f = 0.50),
  adjustcolor("#579ad1", alpha.f = 0.20)
)

sample_colors <- c(
  "Sample" = "#579ad1",
  "Blank"  = "#df6765"
)

#Plot
pdf(
  file.path(
    plot_path,
    "Retention_time_alignment2_positivemap2.pdf"
  ),
  width = 14,
  height = 10
)
par(mfrow = c(2, 1))

plot(
  bpc_raw,
  peakType = "none",
  col = color_mapping,
  main = "Base Peak Chromatogram before Retention Time Alignment"
)
grid()
legend(
  "topright",
  legend = names(sample_colors),
  col = sample_colors,
  pch = 15,
  cex = 0.8,
  title = "Sample Type",
  bg = "white"
)

plot(
  bpc_aligned,
  peakType = "none",
  col = color_mapping,
  main = "Base Peak Chromatogram after Retention Time Alignment"
)
grid()

legend(
  "topright",
  legend = names(sample_colors),
  col = sample_colors,
  pch = 15,
  cex = 0.8,
  title = "Sample Type",
  bg = "white"
)

dev.off()

#####2.g) Correspondence ####

load(
  file.path(
    env_path,
    "myEnvironment3_RT_alignement_positive1.RData"
  )
)
register(SerialParam())

mse1 <- groupChromPeaks(
  mse1,
  param = PeakDensityParam(
    sampleGroups = sampleData(mse1)$ESSENTIAL_sample_id,  
    minFraction = 1/3,
    binSize = 0.01,
    bw = 1.8
  )
)

#to save the environment (useful if very long process):
save.image(
  file = file.path(
    env_path,
    "myEnvironment4_peak_grouping13_positive1.RData"
  )
)
#Choline
eic_choline <- chromatogram(mse1, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70), aggregationFun = "max")
#Adenosine
eic_Adenosine <- chromatogram(mse1, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60), aggregationFun = "max")
#Phosphocholine
eic_Phosphocholine <- chromatogram(mse1, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60), aggregationFun = "max")
#Nicotinic_acid
eic_Nicotinic_acid <- chromatogram(mse1, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60), aggregationFun = "max")
#Arginine
eic_Arginine <- chromatogram(mse1, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60), aggregationFun = "max")
#Val-Leu
eic_Val_Leu <- chromatogram(mse1, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60), aggregationFun = "max")
#Homarine
eic_Homarine <- chromatogram(mse1, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60), aggregationFun = "max")
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(mse1, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60), aggregationFun = "max")
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(mse1, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60), aggregationFun = "max")
#Diphenylamine
eic_Diphenylamine <- chromatogram(mse1, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60), aggregationFun = "max")
#Thr_Leu
eic_Thr_Leu <- chromatogram(mse1, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60), aggregationFun = "max")

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(mse1, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60), aggregationFun = "max")
#Betain_lipid
eic_Betain_lipid <- chromatogram(mse1, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60), aggregationFun = "max")

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(mse1, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60), aggregationFun = "max")

#Tyrosine
eic_Tyrosine <- chromatogram(mse1, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60), aggregationFun = "max")

#Caffeine
eic_Caffeine <- chromatogram(mse1, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60), aggregationFun = "max")
#Vanillin
eic_Vanillin <- chromatogram(mse1, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60), aggregationFun = "max")
#Stachydrine
eic_Stachydrine <- chromatogram(mse1, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60), aggregationFun = "max")

#homocysteine
eic_homocysteine <- chromatogram(mse1, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-60, 60), aggregationFun = "max")
#Guanine
eic_Guanine <- chromatogram(mse1, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60), aggregationFun = "max")
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(mse1, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60), aggregationFun = "max")
#Riboflavin
eic_Riboflavin <- chromatogram(mse1, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60), aggregationFun = "max")
#Cholesterol
eic_Cholesterol <- chromatogram(mse1, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60), aggregationFun = "max")
#Ile_Lys
eic_Ile_Lys <- chromatogram(mse1, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60), aggregationFun = "max")

#Loliolide
eic_Loliolide <- chromatogram(mse1, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60), aggregationFun = "max")
#myrcene
eic_Myrcene <- chromatogram(mse1, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60), aggregationFun = "max")
#Astaxanthin
eic_Astaxanthin <- chromatogram(mse1, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60), aggregationFun = "max")
#LUMICHROME
eic_Lumichrome <- chromatogram(mse1, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60), aggregationFun = "max")
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(mse1, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60), aggregationFun = "max")
#BETAINE
eic_Betain <- chromatogram(mse1, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60), aggregationFun = "max")

#BiakamideC,D
eic_Biakamide <- chromatogram(mse1, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60), aggregationFun = "max")
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(mse1, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60), aggregationFun = "max")
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(mse1, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60), aggregationFun = "max")
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(mse1, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60), aggregationFun = "max")
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(mse1, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60), aggregationFun = "max")

#Shinorine
eic_Shinorine  <- chromatogram(mse1, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60), aggregationFun = "max")
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(mse1, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60), aggregationFun = "max")
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(mse1, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60), aggregationFun = "max")
#Usujirene
eic_Usujirene <- chromatogram(mse1, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60), aggregationFun = "max")
#ANTHRANILATE
eic_Anthranilate <- chromatogram(mse1, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60), aggregationFun = "max")
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(mse1, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#Ala-Phe
eic_Ala_Phe <- chromatogram(mse1, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60), aggregationFun = "max")
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(mse1, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60), aggregationFun = "max")
#Octinoxate
eic_Octinoxate <- chromatogram(mse1, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60), aggregationFun = "max")
#Val-Trp
eic_Val_Trp <- chromatogram(mse1, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60), aggregationFun = "max")
#Inosine
eic_Inosine <- chromatogram(mse1, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60), aggregationFun = "max")
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(mse1, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60), aggregationFun = "max")
#Smenamide A
eic_Smenamide <- chromatogram(mse1, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60), aggregationFun = "max")
#Ser-Ile
eic_Ser_Ile <- chromatogram(mse1, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60), aggregationFun = "max")
#Phe-Glu
eic_Phe_Glu <- chromatogram(mse1, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60), aggregationFun = "max")
#Petasol
eic_Petasol <- chromatogram(mse1, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60), aggregationFun = "max")
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(mse1, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(mse1, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60), aggregationFun = "max")
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(mse1, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60), aggregationFun = "max")
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(mse1, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60), aggregationFun = "max")


save.image(
  file = file.path(
    env_path,
    "myEnvironment4_peak_grouping_EIC13_positive1.RData"
  )
)


#PDF
pdf(
  file.path(
    plot_path,
    "Peak_grouping13_positive3.pdf"
  ),
  width = 14,
  height = 10
)

par(mfrow = c(4, 5))
peak_bg_color <- adjustcolor("cyan4", alpha.f = 0.3)

eic_list <- list(
  eic_choline = expression(paste(italic(m/z), " 104.11")),
  eic_Adenosine = expression(paste(italic(m/z), " 268.10")),
  eic_Phosphocholine = expression(paste(italic(m/z), " 184.07")),
  eic_Nicotinic_acid = expression(paste(italic(m/z), " 124.04")),
  eic_Arginine = expression(paste(italic(m/z), " 175.12")),
  eic_Val_Leu = expression(paste(italic(m/z), " 231.17")),
  eic_Homarine = expression(paste(italic(m/z), " 138.06")),
  eic_Salicylic_acid = expression(paste(italic(m/z), " 121.03")),
  eic_Lauryl_diethanolamide = expression(paste(italic(m/z), " 288.25")),
  eic_Diphenylamine = expression(paste(italic(m/z), " 170.10")),
  eic_Thr_Leu = expression(paste(italic(m/z), " 233.15")),
  eic_CocamidoprpylBetaine = expression(paste(italic(m/z), " 343.30")),
  eic_Betain_lipid = expression(paste(italic(m/z), " 490.37")),
  eic_p_coumaric_acid = expression(paste(italic(m/z), " 147.04")),
  eic_Tyrosine = expression(paste(italic(m/z), " 182.08")),
  eic_Caffeine = expression(paste(italic(m/z), " 195.09")),
  eic_Vanillin = expression(paste(italic(m/z), " 153.06")),
  eic_Stachydrine = expression(paste(italic(m/z), " 144.10")),
  eic_homocysteine = expression(paste(italic(m/z), " 385.13")),
  eic_Guanine = expression(paste(italic(m/z), " 152.06")),
  eic_Dinoxanthin = expression(paste(italic(m/z), " 641.42")),
  eic_Riboflavin = expression(paste(italic(m/z), " 377.15")),
  eic_Cholesterol = expression(paste(italic(m/z), " 369.35")),
  eic_Ile_Lys = expression(paste(italic(m/z), " 260.20")),
  eic_Loliolide = expression(paste(italic(m/z), " 197.12")),
  eic_Myrcene = expression(paste(italic(m/z), " 159.12")),
  eic_Astaxanthin = expression(paste(italic(m/z), " 597.39")),
  eic_Lumichrome = expression(paste(italic(m/z), " 243.09")),
  eic_Diadinoxanthin = expression(paste(italic(m/z), " 581.40")),
  eic_Betain = expression(paste(italic(m/z), " 118.09")),
  eic_Biakamide = expression(paste(italic(m/z), " 524.27")), #Biakamide
  eic_Fucoxanthin = expression(paste(italic(m/z), " 659.43")),
  eic_Arachidonic_acid = expression(paste(italic(m/z), " 305.25")),
  eic_Phaeophorbide = expression(paste(italic(m/z), " 593.27")),
  eic_Crassostreaxanthin = expression(paste(italic(m/z), " 599.41")),
  eic_Shinorine = expression(paste(italic(m/z), " 333.13")),
  eic_Pinolenic_acid = expression(paste(italic(m/z), " 161.10")),
  eic_Trimethyllysine = expression(paste(italic(m/z), " 189.16")),
  eic_Usujirene = expression(paste(italic(m/z), " 285.14")),
  eic_Anthranilate = expression(paste(italic(m/z), " 120.04")),
  eic_Succinoadenosine = expression(paste(italic(m/z), " 384.12")),
  eic_Ala_Phe = expression(paste(italic(m/z), " 237.12")),
  eic_Diphenylphosphate = expression(paste(italic(m/z), " 251.05")),
  eic_Octinoxate = expression(paste(italic(m/z), " 291.19")),
  eic_Val_Trp = expression(paste(italic(m/z), " 304.17")),
  eic_Inosine = expression(paste(italic(m/z), " 269.09")),
  eic_Methyl_coumarin = expression(paste(italic(m/z), " 161.06")),
  eic_Smenamide = expression(paste(italic(m/z), " 501.25")),
  eic_Ser_Ile = expression(paste(italic(m/z), " 219.13")),
  eic_Phe_Glu = expression(paste(italic(m/z), " 295.13")),
  eic_Petasol = expression(paste(italic(m/z), " 235.17")),
  eic_Butyrylcarnitine = expression(paste(italic(m/z), " 232.15")),
  eic_palmitoyl_carnitine = expression(paste(italic(m/z), " 400.34")),
  eic_Benzothiazolesulfonic_acid = expression(paste(italic(m/z), " 215.98")),
  eic_Lauroylcarnitine = expression(paste(italic(m/z), " 344.28"))
)

for (eic_name in names(eic_list)) {
  tryCatch({
    plotChromPeakDensity(get(eic_name), main = eic_list[[eic_name]], peakBg = peak_bg_color, simulate = FALSE)
  }, error = function(e) {
    message(sprintf("Error with %s: %s", eic_name, e$message))
  })
}

dev.off()

#####2.h) Gap filling ####

#Fill in the missing values in the whole dataset
mse1 <- fillChromPeaks(mse1, param = ChromPeakAreaParam(), chunkSize = 5)

save.image(
  file = file.path(
    env_path,
    "myEnvironment5_gap_filling_positive1.RData"
  )
)

#Choline
eic_choline <- chromatogram(mse1, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70), aggregationFun = "max")
#Adenosine
eic_Adenosine <- chromatogram(mse1, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60), aggregationFun = "max")
#Phosphocholine
eic_Phosphocholine <- chromatogram(mse1, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60), aggregationFun = "max")
#Nicotinic_acid
eic_Nicotinic_acid <- chromatogram(mse1, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60), aggregationFun = "max")
#Arginine
eic_Arginine <- chromatogram(mse1, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60), aggregationFun = "max")
#Val-Leu
eic_Val_Leu <- chromatogram(mse1, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60), aggregationFun = "max")
#Homarine
eic_Homarine <- chromatogram(mse1, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60), aggregationFun = "max")
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(mse1, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60), aggregationFun = "max")
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(mse1, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60), aggregationFun = "max")
#Diphenylamine
eic_Diphenylamine <- chromatogram(mse1, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60), aggregationFun = "max")
#Thr_Leu
eic_Thr_Leu <- chromatogram(mse1, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60), aggregationFun = "max")

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(mse1, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60), aggregationFun = "max")
#Betain_lipid
eic_Betain_lipid <- chromatogram(mse1, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60), aggregationFun = "max")

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(mse1, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60), aggregationFun = "max")

#Tyrosine
eic_Tyrosine <- chromatogram(mse1, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60), aggregationFun = "max")

#Caffeine
eic_Caffeine <- chromatogram(mse1, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60), aggregationFun = "max")
#Vanillin
eic_Vanillin <- chromatogram(mse1, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60), aggregationFun = "max")
#Stachydrine
eic_Stachydrine <- chromatogram(mse1, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60), aggregationFun = "max")

#homocysteine
eic_homocysteine <- chromatogram(mse1, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-60, 60), aggregationFun = "max")
#Guanine
eic_Guanine <- chromatogram(mse1, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60), aggregationFun = "max")
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(mse1, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60), aggregationFun = "max")
#Riboflavin
eic_Riboflavin <- chromatogram(mse1, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60), aggregationFun = "max")
#Cholesterol
eic_Cholesterol <- chromatogram(mse1, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60), aggregationFun = "max")
#Ile_Lys
eic_Ile_Lys <- chromatogram(mse1, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60), aggregationFun = "max")

#Loliolide
eic_Loliolide <- chromatogram(mse1, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60), aggregationFun = "max")
#myrcene
eic_Myrcene <- chromatogram(mse1, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60), aggregationFun = "max")
#Astaxanthin
eic_Astaxanthin <- chromatogram(mse1, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60), aggregationFun = "max")
#LUMICHROME
eic_Lumichrome <- chromatogram(mse1, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60), aggregationFun = "max")
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(mse1, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60), aggregationFun = "max")
#BETAINE
eic_Betain <- chromatogram(mse1, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60), aggregationFun = "max")

#BiakamideC,D
eic_Biakamide <- chromatogram(mse1, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60), aggregationFun = "max")
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(mse1, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60), aggregationFun = "max")
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(mse1, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60), aggregationFun = "max")
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(mse1, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60), aggregationFun = "max")
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(mse1, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60), aggregationFun = "max")

#Shinorine
eic_Shinorine  <- chromatogram(mse1, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60), aggregationFun = "max")
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(mse1, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60), aggregationFun = "max")
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(mse1, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60), aggregationFun = "max")
#Usujirene
eic_Usujirene <- chromatogram(mse1, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60), aggregationFun = "max")
#ANTHRANILATE
eic_Anthranilate <- chromatogram(mse1, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60), aggregationFun = "max")
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(mse1, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#Ala-Phe
eic_Ala_Phe <- chromatogram(mse1, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60), aggregationFun = "max")
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(mse1, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60), aggregationFun = "max")
#Octinoxate
eic_Octinoxate <- chromatogram(mse1, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60), aggregationFun = "max")
#Val-Trp
eic_Val_Trp <- chromatogram(mse1, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60), aggregationFun = "max")
#Inosine
eic_Inosine <- chromatogram(mse1, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60), aggregationFun = "max")
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(mse1, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60), aggregationFun = "max")
#Smenamide A
eic_Smenamide <- chromatogram(mse1, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60), aggregationFun = "max")
#Ser-Ile
eic_Ser_Ile <- chromatogram(mse1, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60), aggregationFun = "max")
#Phe-Glu
eic_Phe_Glu <- chromatogram(mse1, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60), aggregationFun = "max")
#Petasol
eic_Petasol <- chromatogram(mse1, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60), aggregationFun = "max")
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(mse1, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(mse1, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60), aggregationFun = "max")
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(mse1, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60), aggregationFun = "max")
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(mse1, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60), aggregationFun = "max")
save.image(
  file = file.path(
    env_path,
    "myEnvironment4_gap_fillingg_EIC13_positive1.RData"
  )
)

#PDF
pdf(
  file.path(
    plot_path,
    "Peak_gap_filling_positive1.pdf"
  ),
  width = 14,
  height = 10
)

par(mfrow = c(4, 5))
peak_bg_color <- adjustcolor("cyan4", alpha.f = 0.3)
eic_list <- list(
  eic_choline = expression(paste(italic(m/z), " 104.11")),
  eic_Adenosine = expression(paste(italic(m/z), " 268.10")),
  eic_Phosphocholine = expression(paste(italic(m/z), " 184.07")),
  eic_Nicotinic_acid = expression(paste(italic(m/z), " 124.04")),
  eic_Arginine = expression(paste(italic(m/z), " 175.12")),
  eic_Val_Leu = expression(paste(italic(m/z), " 231.17")),
  eic_Homarine = expression(paste(italic(m/z), " 138.06")),
  eic_Salicylic_acid = expression(paste(italic(m/z), " 121.03")),
  eic_Lauryl_diethanolamide = expression(paste(italic(m/z), " 288.25")),
  eic_Diphenylamine = expression(paste(italic(m/z), " 170.10")),
  eic_Thr_Leu = expression(paste(italic(m/z), " 233.15")),
  eic_CocamidoprpylBetaine = expression(paste(italic(m/z), " 343.30")),
  eic_Betain_lipid = expression(paste(italic(m/z), " 490.37")),
  eic_p_coumaric_acid = expression(paste(italic(m/z), " 147.04")),
  eic_Tyrosine = expression(paste(italic(m/z), " 182.08")),
  eic_Caffeine = expression(paste(italic(m/z), " 195.09")),
  eic_Vanillin = expression(paste(italic(m/z), " 153.06")),
  eic_Stachydrine = expression(paste(italic(m/z), " 144.10")),
  eic_homocysteine = expression(paste(italic(m/z), " 385.13")),
  eic_Guanine = expression(paste(italic(m/z), " 152.06")),
  eic_Dinoxanthin = expression(paste(italic(m/z), " 641.42")),
  eic_Riboflavin = expression(paste(italic(m/z), " 377.15")),
  eic_Cholesterol = expression(paste(italic(m/z), " 369.35")),
  eic_Ile_Lys = expression(paste(italic(m/z), " 260.20")),
  eic_Loliolide = expression(paste(italic(m/z), " 197.12")),
  eic_Myrcene = expression(paste(italic(m/z), " 159.12")),
  eic_Astaxanthin = expression(paste(italic(m/z), " 597.39")),
  eic_Lumichrome = expression(paste(italic(m/z), " 243.09")),
  eic_Diadinoxanthin = expression(paste(italic(m/z), " 581.40")),
  eic_Betain = expression(paste(italic(m/z), " 118.09")),
  eic_Biakamide = expression(paste(italic(m/z), " 524.27")), #Biakamide
  eic_Fucoxanthin = expression(paste(italic(m/z), " 659.43")),
  eic_Arachidonic_acid = expression(paste(italic(m/z), " 305.25")),
  eic_Phaeophorbide = expression(paste(italic(m/z), " 593.27")),
  eic_Crassostreaxanthin = expression(paste(italic(m/z), " 599.41")),
  eic_Shinorine = expression(paste(italic(m/z), " 333.13")),
  eic_Pinolenic_acid = expression(paste(italic(m/z), " 161.10")),
  eic_Trimethyllysine = expression(paste(italic(m/z), " 189.16")),
  eic_Usujirene = expression(paste(italic(m/z), " 285.14")),
  eic_Anthranilate = expression(paste(italic(m/z), " 120.04")),
  eic_Succinoadenosine = expression(paste(italic(m/z), " 384.12")),
  eic_Ala_Phe = expression(paste(italic(m/z), " 237.12")),
  eic_Diphenylphosphate = expression(paste(italic(m/z), " 251.05")),
  eic_Octinoxate = expression(paste(italic(m/z), " 291.19")),
  eic_Val_Trp = expression(paste(italic(m/z), " 304.17")),
  eic_Inosine = expression(paste(italic(m/z), " 269.09")),
  eic_Methyl_coumarin = expression(paste(italic(m/z), " 161.06")),
  eic_Smenamide = expression(paste(italic(m/z), " 501.25")),
  eic_Ser_Ile = expression(paste(italic(m/z), " 219.13")),
  eic_Phe_Glu = expression(paste(italic(m/z), " 295.13")),
  eic_Petasol = expression(paste(italic(m/z), " 235.17")),
  eic_Butyrylcarnitine = expression(paste(italic(m/z), " 232.15")),
  eic_palmitoyl_carnitine = expression(paste(italic(m/z), " 400.34")),
  eic_Benzothiazolesulfonic_acid = expression(paste(italic(m/z), " 215.98")),
  eic_Lauroylcarnitine = expression(paste(italic(m/z), " 344.28"))
)

for (eic_name in names(eic_list)) {
  tryCatch({
    plotChromPeakDensity(get(eic_name), main = eic_list[[eic_name]], peakBg = peak_bg_color, simulate = FALSE)
  }, error = function(e) {
    message(sprintf("Error with %s: %s", eic_name, e$message))
  })
}


dev.off()


#####2.i) Save the data ####

#Alabaster (Structured archival using HDF5 and JSON):
#Compatible with MsExperiment, Spectra, XcmsExperiment, etc.
saveMsObject(
  mse1,
  AlabasterParam(
    path = file.path(
      out_path,
      "Alabaster_preprocessed_MS1pos"
    )
  )
)
#Plain Text (Tab-delimited export/import for key objects) :
#MsBackendMzR, Spectra (from Spectra)
#MsExperiment (from MsExperiment)
#XcmsExperiment (from xcms)
saveMsObject(
  mse1,
  PlainTextParam(
    path = file.path(
      out_path,
      "PlainText_preprocessed_MS1pos"
    )
  )
)
#mzTab-M Export (HUPO PSI metabolomics standard):
param <- mzTabParam(
  studyId = "MzTab_preprocessed_MS1pos",
  polarity = "positive",
  sampleDataColumn = colnames(sampleData(mse1)),
  path = out_path
)

saveMsObject(mse1, param)


#Peak table
#Export peak area quantification table. 
featuresDef <- featureDefinitions(mse1)
featuresIntensities <- featureValues(mse1, value = "into")
dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
head(dataTable)
write.table(
  dataTable,
  "feature_table_MS1_positive1.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


