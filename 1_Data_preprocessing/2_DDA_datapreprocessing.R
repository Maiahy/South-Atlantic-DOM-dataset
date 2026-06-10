#XCMS code TARA for the DRACO cluster

#this code is based on the work of different groups, for more information please consult:
#https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html
#https://github.com/DorresteinLaboratory/XCMS3_FeatureBasedMN/blob/master/XCMS3_Preprocessing.Rmd
#https://pietrofranceschi.github.io/Metabolomics_lectures/
#https://rformassspectrometry.github.io/Metabonaut/articles/a-end-to-end-untargeted-metabolomics.html

# Define main path
library(here)
base_path <- here()

# Define sub-paths
raw_data_path <- file.path(base_path, "1_Raw_data")
metadata_path <- file.path(base_path, "Metadata_final_with_operators.txt")
out_path <- file.path(base_path, "3_Output_files")
env_path <- file.path(out_path, "1_Environments")
plot_path <- file.path(out_path, "2_Quality_plots")

# Create directories (you need to have 1_Raw_data and Metadata_final_with_operators.txt created)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(env_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# Set personal library path
#.libPaths("/home/re73voq/R/x86_64-redhat-linux-gnu-library/4.2")

# Load required packages
source(file.path(base_path, "customFunctions.R")) # https://github.com/jorainer/xcms-gnps-tools/blob/master/customFunctions.R , accessed on the 07/11/2024
library(xcms)
library(CAMERA)
library(RColorBrewer)
library(magrittr)
library(SummarizedExperiment)
library(MSnbase)
library(janitor)
library(dplyr)
library(BiocParallel)


#++++++++++++++++++++++++++++++++++++++++++++++++++++
##### 0) Import data (readMSData)####
#++++++++++++++++++++++++++++++++++++++++++++++++++++

# Register parallel backend
register(MulticoreParam(45))

# Read metadata 
metadata1 <- read.table(metadata_path, sep = "\t") %>%
  row_to_names(row_number = 1) 

# List files in the raw data directory
mydata1 <- list.files(path = raw_data_path, full.names = TRUE) 

#Keep only the files from DDA dataset
inclusion_patterns <- c("_B1_", "_B2_", "_B3_", "_B4_", "_B5_")
inclusion_patterns <- paste(inclusion_patterns, collapse = "|")
filtered_files <- mydata1[grepl(inclusion_patterns, mydata1, ignore.case = TRUE)]

#Create your MsExperiment object
pd <- data.frame(file_name = basename(filtered_files))

# Read metadata and join with file list
pd <- left_join(pd, metadata1, by = "file_name")

#we need to change once more the directory in order to make the files
#accessibles
setwd(base_path)
raw_data <- readMSData(filtered_files, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk") #create my raw data object, could be long
raw_data


#Choose polarity mode to work with: 1 for positive, 0 for negative
polarity_mode <- 1  

#Filter raw data to keep only the selected polarity
raw_data <- filterPolarity(raw_data, polarity = polarity_mode)

#Create the table
ms_polarity_table <- table(msLevel(raw_data), polarity(raw_data))

#Convert it into a data frame 
ms_polarity_df <- as.data.frame(ms_polarity_table)

#Save it as a .csv file
write.csv(ms_polarity_df, file = file.path(out_path, paste0("DDAmsLevel_vs_polarity_postive1", ".csv")), 
          row.names = FALSE)

#to save the environment (useful if very long process):
save.image(
  file = file.path(
    env_path,
    "myEnvironment1_Data_Loading1.RData"
  )
)

##### 1) Peak picking####

#Set the parameters
#noise = 2E5 estimated with mzmine3 (v3.2.8) with the pooled QC n°25
#peakwidth optimised according the same method
cwp <- CentWaveParam(snthresh = 5, noise = 2E5, peakwidth = c(2, 20), 
                     ppm = 15, prefilter = c(3, 3000), integrate = 2)

#Perform the chromatographic peak detection using *centWave*.
processedData <- findChromPeaks(raw_data, param = cwp)

chromPeaks(processedData) |>
  head()

#Check for some compounds of interest:

#Choline
eic_choline <- chromatogram(processedData, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70))
#Adenosine
eic_Adenosine <- chromatogram(processedData, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-30, 30))
#Phosphocholine
eic_Phosphocholine <- chromatogram(processedData, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-30, 30))
#Nicotinic_acid
eic_Nicotinic_acid <- chromatogram(processedData, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-30, 30))
#Arginine
eic_Arginine <- chromatogram(processedData, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-30, 30))
#Val-Leu
eic_Val_Leu <- chromatogram(processedData, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-30, 30))
#Homarine
eic_Homarine <- chromatogram(processedData, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-30, 30))
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(processedData, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-30, 30))
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(processedData, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-30, 30))
#Diphenylamine
eic_Diphenylamine <- chromatogram(processedData, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-30, 30))
#Thr_Leu
eic_Thr_Leu <- chromatogram(processedData, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-30, 30))

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(processedData, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-30, 30))
#Betain_lipid
eic_Betain_lipid <- chromatogram(processedData, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-30, 30))

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(processedData, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-30, 30))

#Tyrosine
eic_Tyrosine <- chromatogram(processedData, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-30, 30))

#Caffeine
eic_Caffeine <- chromatogram(processedData, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-30, 30))
#Vanillin
eic_Vanillin <- chromatogram(processedData, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-30, 30))
#Stachydrine
eic_Stachydrine <- chromatogram(processedData, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-30, 30))

#homocysteine
eic_homocysteine <- chromatogram(processedData, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-30, 30))
#Guanine
eic_Guanine <- chromatogram(processedData, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-30, 30))
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(processedData, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-30, 30))
#Riboflavin
eic_Riboflavin <- chromatogram(processedData, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-30, 30))
#Cholesterol
eic_Cholesterol <- chromatogram(processedData, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-30, 30))
#Ile_Lys
eic_Ile_Lys <- chromatogram(processedData, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-30, 30))

#Loliolide
eic_Loliolide <- chromatogram(processedData, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-30, 30))
#myrcene
eic_Myrcene <- chromatogram(processedData, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-30, 30))
#Astaxanthin
eic_Astaxanthin <- chromatogram(processedData, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-30, 30))
#LUMICHROME
eic_Lumichrome <- chromatogram(processedData, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-30, 30))
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(processedData, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-30, 30))
#BETAINE
eic_Betain <- chromatogram(processedData, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-30, 30))

#BiakamideC,D
eic_Biakamide <- chromatogram(processedData, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-30, 30))
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(processedData, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-30, 30))
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(processedData, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-30, 30))
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(processedData, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-30, 30))
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(processedData, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-30, 30))

#Shinorine
eic_Shinorine  <- chromatogram(processedData, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-30, 30))
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(processedData, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-30, 30))
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(processedData, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-30, 30))
#Usujirene
eic_Usujirene <- chromatogram(processedData, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-30, 30))
#ANTHRANILATE
eic_Anthranilate <- chromatogram(processedData, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-30, 30))
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(processedData, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-30, 30))
#Ala-Phe
eic_Ala_Phe <- chromatogram(processedData, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-30, 30))
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(processedData, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-30, 30))
#Octinoxate
eic_Octinoxate <- chromatogram(processedData, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-30, 30))
#Val-Trp
eic_Val_Trp <- chromatogram(processedData, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-30, 30))
#Inosine
eic_Inosine <- chromatogram(processedData, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-30, 30))
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(processedData, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-30, 30))
#Smenamide A
eic_Smenamide <- chromatogram(processedData, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-30, 30))
#Ser-Ile
eic_Ser_Ile <- chromatogram(processedData, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-30, 30))
#Phe-Glu
eic_Phe_Glu <- chromatogram(processedData, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-30, 30))
#Petasol
eic_Petasol <- chromatogram(processedData, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-30, 30))
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(processedData, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-30, 30))
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(processedData, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-30, 30))
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(processedData, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-30, 30))
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(processedData, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-30, 30))

library(grDevices)
pdf(
  file.path(
    plot_path,
    "DDAPeak_picking_quality_positive1.pdf"
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

save.image(
  file = file.path(
    env_path,
    "DDAmyEnvironment2_Peak_picking_positive1.RData"
  )
)

#####2) Peak refinement ####

register(SerialParam()) #In case you can't for some reasons

param <- MergeNeighboringPeaksParam(expandRt = 2.5, expandMz = 0.0015,
                                    minProp = 0.75)
processedData <- refineChromPeaks(processedData, param = param)

save.image(
  file = file.path(
    env_path,
    "DDAmyEnvironment2_Peak_refinement_positive1.RData"
  )
)

register(SerialParam()) 
#Choline?
eic_choline <- chromatogram(processedData, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70))
#Adenosine?
eic_Adenosine <- chromatogram(processedData, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60))
#Phosphocholine?
eic_Phosphocholine <- chromatogram(processedData, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60))
#Nicotinic_acid?
eic_Nicotinic_acid <- chromatogram(processedData, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60))
#Arginine
eic_Arginine <- chromatogram(processedData, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60))
#Val-Leu
eic_Val_Leu <- chromatogram(processedData, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60))
#Homarine
eic_Homarine <- chromatogram(processedData, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60))
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(processedData, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60))
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(processedData, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60))
#Diphenylamine
eic_Diphenylamine <- chromatogram(processedData, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60))
#Thr_Leu
eic_Thr_Leu <- chromatogram(processedData, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60))

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(processedData, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60))
#Betain_lipid
eic_Betain_lipid <- chromatogram(processedData, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60))

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(processedData, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60))

#Tyrosine
eic_Tyrosine <- chromatogram(processedData, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60))

#Caffeine
eic_Caffeine <- chromatogram(processedData, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60))
#Vanillin
eic_Vanillin <- chromatogram(processedData, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60))
#Stachydrine
eic_Stachydrine <- chromatogram(processedData, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60))
#homocysteine
eic_homocysteine <- chromatogram(processedData, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-60, 60))
#Guanine
eic_Guanine <- chromatogram(processedData, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60))
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(processedData, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60))
#Riboflavin
eic_Riboflavin <- chromatogram(processedData, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60))
#Cholesterol
eic_Cholesterol <- chromatogram(processedData, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60))
#Ile_Lys
eic_Ile_Lys <- chromatogram(processedData, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60))
#Loliolide
eic_Loliolide <- chromatogram(processedData, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60))
#myrcene
eic_Myrcene <- chromatogram(processedData, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60))
#Astaxanthin
eic_Astaxanthin <- chromatogram(processedData, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60))
#LUMICHROME
eic_Lumichrome <- chromatogram(processedData, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60))
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(processedData, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60))
#BETAINE
eic_Betain <- chromatogram(processedData, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60))
#BiakamideC,D
eic_Biakamide <- chromatogram(processedData, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60))
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(processedData, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60))
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(processedData, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60))
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(processedData, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60))
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(processedData, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60))
#Shinorine
eic_Shinorine  <- chromatogram(processedData, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60))
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(processedData, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60))
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(processedData, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60))
#Usujirene
eic_Usujirene <- chromatogram(processedData, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60))
#ANTHRANILATE
eic_Anthranilate <- chromatogram(processedData, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60))
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(processedData, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#Ala-Phe
eic_Ala_Phe <- chromatogram(processedData, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60))
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(processedData, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60))
#Octinoxate
eic_Octinoxate <- chromatogram(processedData, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60))
#Val-Trp
eic_Val_Trp <- chromatogram(processedData, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60))
#Inosine
eic_Inosine <- chromatogram(processedData, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60))
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(processedData, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60))
#Smenamide A
eic_Smenamide <- chromatogram(processedData, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60))
#Ser-Ile
eic_Ser_Ile <- chromatogram(processedData, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60))
#Phe-Glu
eic_Phe_Glu <- chromatogram(processedData, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60))
#Petasol
eic_Petasol <- chromatogram(processedData, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60))
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(processedData, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60))
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(processedData, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60))
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(processedData, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60))
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(processedData, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60))

save.image(
  file = file.path(
    env_path,
    "DDAmyEnvironment2_Peak_refinement_EIC_positive1.RData"
  )
)


library(grDevices)
pdf(
  file.path(
    plot_path,
    "DDAPeak_refinement_quality_positive2.pdf"
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



### 3) Retention time alignment (adjustRtime). ####

load(
  file.path(
    env_path,
    "DDAmyEnvironment2_Peak_refinement_positive1.RData"
  )
)

register(SerialParam())


processedData1 <- processedData  # Backup before alignment

metadata_check <- data.frame(
  file_name = pData(processedData)$file_name,
  ESSENTIAL_sample_id = pData(processedData)$ESSENTIAL_sample_id,
  ESSENTIAL_sample_type = pData(processedData)$ESSENTIAL_sample_type,
  stringsAsFactors = FALSE
)

#Save to output folder for manual inspection
write.csv(metadata_check,
          file = file.path(out_path, "DDAmetadata_check_sampleType_positive1.csv"),
          row.names = FALSE, quote = TRUE)

#Extract and sort sample metadata by injection order
sample_data <- pData(processedData)
sample_data$Injection_order <- as.numeric(sample_data$Injection_order)
sample_data <- sample_data[order(sample_data$Injection_order), ]
pData(processedData) <- sample_data  # Reassign sorted metadata

#Check injection order
print(pData(processedData)$Injection_order)

#Create a table to save
injection_order_df <- data.frame(
  file_name = sample_data$file_name,
  injection_order = sample_data$Injection_order,
  ESSENTIAL_sample_type = sample_data$ESSENTIAL_sample_type
)

#Save table
write.csv(injection_order_df, file = file.path(out_path, paste0("DDAInjection_order_check_positive1", ".csv")), 
          row.names = FALSE)

processedData <- groupChromPeaks(
  processedData, 
  param = PeakDensityParam(
    sampleGroups = pData(processedData)$ESSENTIAL_sample_id,
    minFraction = 2/3,           
    binSize = 0.02,
    ppm = 10,
    bw = 3))

processedData <- adjustRtime(processedData, 
                             param = PeakGroupsParam(
                               minFraction = 0.8,
                               subset = which(pData(processedData)$ESSENTIAL_sample_id == "InterbatchQC"), 
                               subsetAdjust = "average",
                               span = 0.6))

#Save environment
save.image(
  file = file.path(
    env_path,
    "DDAmyEnvironment3_RT_alignement_positive1.RData"
  )
)
#Quality control figures:

#1)RT adjustment:
sample_colors <- c("QC" = "orangered3", "DDA" = "cyan4")
color_mapping <- sample_colors[processedData$ESSENTIAL_sample_type]
pdf(file.path(plot_path,
              "DDARetention_time_alignment1_positive1_simplified.pdf"),
    width = 14, height = 10)

plotAdjustedRtime(processedData, col = color_mapping,
                  main = "Adjusted retention times across samples")

legend("topright", legend = names(sample_colors),
       col = sample_colors, pch = 15, cex = 0.8,
       title = "Sample Type", bg = "white")

dev.off()

#2)BPCs before and after alignment all samples:
bpc_raw <- chromatogram(processedData1, aggregationFun = "max")
bpc_aligned <- chromatogram(processedData, aggregationFun = "max", include = "none")
pdf(file.path(plot_path,
              "DDARetention_time_alignment2_positive1.pdf"),
    width = 14, height = 10)
par(mfrow = c(2, 1))
plot(bpc_raw, peakType = "none", col = color_mapping, main = "Base Peak Chromatogram before Retention Time Alignment")
grid()
legend("topright", legend = names(sample_colors), col = sample_colors, pch = 15, cex = 0.8, title = "Sample Type", bg = "white")
plot(bpc_aligned, col = color_mapping, main = "Base Peak Chromatogram after Retention Time Alignment")
grid()
dev.off()

#2)BPCs before and after alignment QCs:
qc_idx <- which(pData(processedData)$ESSENTIAL_sample_id == "InterbatchQC")
bpc_raw_qc <- chromatogram(processedData1[qc_idx], aggregationFun = "max")
bpc_aligned_qc <- chromatogram(processedData[qc_idx], aggregationFun = "max", include = "none")
qc_colors <- rep("orangered3", length(qc_idx))
pdf(file.path(plot_path,
              "DDARetention_time_alignment_QC_only_positive1.pdf"),
    width = 14, height = 10)
par(mfrow = c(2, 1))
plot(bpc_raw_qc,
     peakType = "none",
     col = qc_colors,
     main = "BPC before RT alignment (study reference material)")
grid()
plot(bpc_aligned_qc,
     col = qc_colors,
     main = "BPC after RT alignment (study reference material)")
grid()
dev.off()
register(SerialParam())

#3)EIC before and after:

compounds <- data.frame(
  name = c("Choline","Adenosine","Phosphocholine","Nicotinic_acid",
           "Arginine","Val_Leu","Homarine","Salicylic_acid",
           "Lauryl_diethanolamide","Diphenylamine","Thr_Leu",
           "CocamidopropylBetaine","Betain_lipid","p_coumaric_acid",
           "Tyrosine","Caffeine","Vanillin","Stachydrine",
           "Homocysteine","Guanine","Dinoxanthin","Riboflavin",
           "Cholesterol","Ile_Lys","Loliolide","Myrcene",
           "Astaxanthin","Lumichrome","Diadinoxanthin","Betaine",
           "Biakamide","Fucoxanthin","Arachidonic_acid",
           "Phaeophorbide","Crassostreaxanthin","Shinorine",
           "Pinolenic_acid","Trimethyllysine","Usujirene",
           "Anthranilate","Succinoadenosine","Ala_Phe",
           "Diphenylphosphate","Octinoxate","Val_Trp",
           "Inosine","Methyl_coumarin","Smenamide",
           "Ser_Ile","Phe_Glu","Petasol","Butyrylcarnitine",
           "Palmitoyl_carnitine","Benzothiazolesulfonic_acid",
           "Lauroylcarnitine"),
  
  mz = c(104.11,268.10,184.07,124.04,
         175.12,231.17,138.06,121.03,
         288.25,170.10,233.15,
         343.30,490.37,147.04,
         182.08,195.09,153.06,144.10,
         385.13,152.06,641.42,377.15,
         369.35,260.20,197.12,159.12,
         597.39,243.09,581.40,118.09,
         524.27,659.43,305.25,
         593.27,599.41,333.13,
         161.10,189.16,285.14,
         120.04,384.12,237.12,
         251.05,291.19,304.17,
         269.09,161.06,501.25,
         219.13,295.13,235.17,232.15,
         400.34,215.98,344.28),
  
  rt = c(60,90.36,397.18,66.45,
         53.35,183.02,60.21,346.38,
         398.23,421.28,170.58,
         340.58,436.71,75.29,
         144.09,199.16,243.39,60.84,
         72.68,110.74,618.52,197.20,
         650.76,151.55,259.28,388.85,
         606.32,256.55,567.38,60.21,
         437.32,526.04,462.80,
         566.14,563.04,62.44,
         277.36,55.49,148.69,
         81.32,168.22,183.02,
         521.71,438.85,209.34,
         92.87,317.59,432.08,
         153.19,161.04,453.49,168.22,
         476.56,211.42,369.72)
)


pdf(file.path(plot_path,
              "DDA_EIC_before_after_alignment.pdf"),
    width = 10, height = 8, onefile = TRUE)

for (i in seq_len(nrow(compounds))) {
  
  mz_range <- compounds$mz[i] + c(-0.01, 0.01)
  rt_range <- compounds$rt[i] + c(-60, 60)
  
  #Avant alignement
  eic_raw <- chromatogram(processedData1,
                          mz = mz_range,
                          rt = rt_range,
                          aggregationFun = "max")
  
  #Après alignement
  eic_aligned <- chromatogram(processedData,
                              mz = mz_range,
                              rt = rt_range,
                              aggregationFun = "max")
  
  #Nouvelle page
  par(mfrow = c(2, 1))
  
  plot(eic_raw,
       peakType = "none",
       col = "grey40",
       main = paste0("m/z ", compounds$mz[i],
                     "  (Before alignment)"))
  grid()
  
  
  plot(eic_aligned,
       peakType = "none",
       col = "steelblue",
       main = paste0("m/z ", compounds$mz[i],
                     "  (After alignment)"))
  grid()
}

dev.off()


##### 4) Correspondence ####

register(SerialParam())

#Define the parameters for the *peak density*-based peak grouping (correspondence
# analysis).

# Define the parameters for peak density-based peak grouping (correspondence analysis)
processedData <- groupChromPeaks(
  processedData,
  param = PeakDensityParam(
    sampleGroups = pData(processedData)$ESSENTIAL_sample_id,  
    minFraction = 1/3,
    binSize = 0.01,
    bw = 1.8
  )
)

#Choline
eic_choline <- chromatogram(processedData, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70), aggregationFun = "max")
#Adenosine
eic_Adenosine <- chromatogram(processedData, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-30, 30), aggregationFun = "max")
#Phosphocholine
eic_Phosphocholine <- chromatogram(processedData, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-30, 30), aggregationFun = "max")
#Nicotinic_acid
eic_Nicotinic_acid <- chromatogram(processedData, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-30, 30), aggregationFun = "max")
#Arginine
eic_Arginine <- chromatogram(processedData, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-30, 30), aggregationFun = "max")
#Val-Leu
eic_Val_Leu <- chromatogram(processedData, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-30, 30), aggregationFun = "max")
#Homarine
eic_Homarine <- chromatogram(processedData, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-30, 30), aggregationFun = "max")
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(processedData, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-30, 30), aggregationFun = "max")
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(processedData, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-30, 30), aggregationFun = "max")
#Diphenylamine
eic_Diphenylamine <- chromatogram(processedData, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-30, 30), aggregationFun = "max")
#Thr_Leu
eic_Thr_Leu <- chromatogram(processedData, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-30, 30), aggregationFun = "max")

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(processedData, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-30, 30), aggregationFun = "max")
#Betain_lipid
eic_Betain_lipid <- chromatogram(processedData, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-30, 30), aggregationFun = "max")

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(processedData, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-30, 30), aggregationFun = "max")

#Tyrosine
eic_Tyrosine <- chromatogram(processedData, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-30, 30), aggregationFun = "max")

#Caffeine
eic_Caffeine <- chromatogram(processedData, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-30, 30), aggregationFun = "max")
#Vanillin
eic_Vanillin <- chromatogram(processedData, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-30, 30), aggregationFun = "max")
#Stachydrine
eic_Stachydrine <- chromatogram(processedData, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-30, 30), aggregationFun = "max")

#homocysteine
eic_homocysteine <- chromatogram(processedData, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-30, 30), aggregationFun = "max")
#Guanine
eic_Guanine <- chromatogram(processedData, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-30, 30), aggregationFun = "max")
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(processedData, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-30, 30), aggregationFun = "max")
#Riboflavin
eic_Riboflavin <- chromatogram(processedData, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-30, 30), aggregationFun = "max")
#Cholesterol
eic_Cholesterol <- chromatogram(processedData, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-30, 30), aggregationFun = "max")
#Ile_Lys
eic_Ile_Lys <- chromatogram(processedData, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-30, 30), aggregationFun = "max")

#Loliolide
eic_Loliolide <- chromatogram(processedData, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-30, 30), aggregationFun = "max")
#myrcene
eic_Myrcene <- chromatogram(processedData, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-30, 30), aggregationFun = "max")
#Astaxanthin
eic_Astaxanthin <- chromatogram(processedData, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-30, 30), aggregationFun = "max")
#LUMICHROME
eic_Lumichrome <- chromatogram(processedData, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-30, 30), aggregationFun = "max")
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(processedData, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-30, 30), aggregationFun = "max")
#BETAINE
eic_Betain <- chromatogram(processedData, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-30, 30), aggregationFun = "max")

#BiakamideC,D
eic_Biakamide <- chromatogram(processedData, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-30, 30), aggregationFun = "max")
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(processedData, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-30, 30), aggregationFun = "max")
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(processedData, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-30, 30), aggregationFun = "max")
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(processedData, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-30, 30), aggregationFun = "max")
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(processedData, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-30, 30), aggregationFun = "max")

#Shinorine
eic_Shinorine  <- chromatogram(processedData, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-30, 30), aggregationFun = "max")
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(processedData, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-30, 30), aggregationFun = "max")
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(processedData, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-30, 30), aggregationFun = "max")
#Usujirene
eic_Usujirene <- chromatogram(processedData, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-30, 30), aggregationFun = "max")
#ANTHRANILATE
eic_Anthranilate <- chromatogram(processedData, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-30, 30), aggregationFun = "max")
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(processedData, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-30, 30), aggregationFun = "max")
#Ala-Phe
eic_Ala_Phe <- chromatogram(processedData, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-30, 30), aggregationFun = "max")
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(processedData, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-30, 30), aggregationFun = "max")
#Octinoxate
eic_Octinoxate <- chromatogram(processedData, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-30, 30), aggregationFun = "max")
#Val-Trp
eic_Val_Trp <- chromatogram(processedData, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-30, 30), aggregationFun = "max")
#Inosine
eic_Inosine <- chromatogram(processedData, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-30, 30), aggregationFun = "max")
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(processedData, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-30, 30), aggregationFun = "max")
#Smenamide A
eic_Smenamide <- chromatogram(processedData, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-30, 30), aggregationFun = "max")
#Ser-Ile
eic_Ser_Ile <- chromatogram(processedData, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-30, 30), aggregationFun = "max")
#Phe-Glu
eic_Phe_Glu <- chromatogram(processedData, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-30, 30), aggregationFun = "max")
#Petasol
eic_Petasol <- chromatogram(processedData, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-30, 30), aggregationFun = "max")
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(processedData, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-30, 30), aggregationFun = "max")
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(processedData, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-30, 30), aggregationFun = "max")
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(processedData, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-30, 30), aggregationFun = "max")
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(processedData, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-30, 30), aggregationFun = "max")

# Open PDF device
pdf(file.path(plot_path,
              "DDAPeak_grouping_positive1.pdf"),
    width = 14, height = 10)

par(mfrow = c(4, 5))
peak_bg_color <- adjustcolor("cyan4", alpha.f = 0.3)

# Create a vector of the EIC objects and their corresponding titles
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


# Loop through each EIC and plot it, using tryCatch to handle errors
for (eic_name in names(eic_list)) {
  tryCatch({
    plotChromPeakDensity(get(eic_name), main = eic_list[[eic_name]], peakBg = peak_bg_color, simulate = FALSE)
  }, error = function(e) {
    message(sprintf("Error with %s: %s", eic_name, e$message))
  })
}

# Close the PDF device
dev.off()

#to save the environment (useful if very long process):
save.image(
  file = file.path(
    env_path,
    "DDAmyEnvironment4_peak_grouping_positive1.RData"
  )
)


#### 5) Gap filling (fillChromPeaks). ####

register(SerialParam())

## fill missing peaks
processed_Data <- fillChromPeaks(processedData, param = ChromPeakAreaParam())
save.image(file = file.path(env_path,
                            "DDAmyEnvironment5_fillChromPeak1.RData"))
#Choline
eic_choline <- chromatogram(processed_Data, mz = 104.11 + c(-0.01, 0.01),
                            rt = c(50, 70), aggregationFun = "max")
#Adenosine
eic_Adenosine <- chromatogram(processed_Data, mz = 268.10 + c(-0.01, 0.01),
                              rt = 90.36 + c(-60, 60), aggregationFun = "max")
#Phosphocholine
eic_Phosphocholine <- chromatogram(processed_Data, mz = 	184.07 + c(-0.01, 0.01),
                                   rt = 397.18 + c(-60, 60), aggregationFun = "max")
#Nicotinic_acid
eic_Nicotinic_acid <- chromatogram(processed_Data, mz = 	124.04 + c(-0.01, 0.01),
                                   rt = 66.45 + c(-60, 60), aggregationFun = "max")
#Arginine
eic_Arginine <- chromatogram(processed_Data, mz = 	175.12 + c(-0.01, 0.01),
                             rt = 53.35 + c(-60, 60), aggregationFun = "max")
#Val-Leu
eic_Val_Leu <- chromatogram(processed_Data, mz = 	231.17 + c(-0.01, 0.01),
                            rt = 183.02 + c(-60, 60), aggregationFun = "max")
#Homarine
eic_Homarine <- chromatogram(processed_Data, mz = 	138.06 + c(-0.01, 0.01),
                             rt = 60.21 + c(-60, 60), aggregationFun = "max")
#Salicylic_acid
eic_Salicylic_acid <- chromatogram(processed_Data, mz = 	121.03 + c(-0.01, 0.01),
                                   rt = 346.38 + c(-60, 60), aggregationFun = "max")
#Lauryl_diethanolamide
eic_Lauryl_diethanolamide <- chromatogram(processed_Data, mz = 	288.25 + c(-0.01, 0.01),
                                          rt = 398.23 + c(-60, 60), aggregationFun = "max")
#Diphenylamine
eic_Diphenylamine <- chromatogram(processed_Data, mz = 	170.10 + c(-0.01, 0.01),
                                  rt = 421.28 + c(-60, 60), aggregationFun = "max")
#Thr_Leu
eic_Thr_Leu <- chromatogram(processed_Data, mz = 	233.15 + c(-0.01, 0.01),
                            rt = 170.58 + c(-60, 60), aggregationFun = "max")

#CocamidoprpylBetaine
eic_CocamidoprpylBetaine <- chromatogram(processed_Data, mz = 	343.30 + c(-0.01, 0.01),
                                         rt = 340.58 + c(-60, 60), aggregationFun = "max")
#Betain_lipid
eic_Betain_lipid <- chromatogram(processed_Data, mz = 	490.37 + c(-0.01, 0.01),
                                 rt = 436.71 + c(-60, 60), aggregationFun = "max")

#p_coumaric_acid
eic_p_coumaric_acid <- chromatogram(processed_Data, mz = 	147.04 + c(-0.01, 0.01),
                                    rt = 75.29 + c(-60, 60), aggregationFun = "max")

#Tyrosine
eic_Tyrosine <- chromatogram(processed_Data, mz = 	182.08 + c(-0.01, 0.01),
                             rt = 144.09 + c(-60, 60), aggregationFun = "max")

#Caffeine
eic_Caffeine <- chromatogram(processed_Data, mz = 	195.09 + c(-0.01, 0.01),
                             rt = 199.16 + c(-60, 60), aggregationFun = "max")
#Vanillin
eic_Vanillin <- chromatogram(processed_Data, mz = 	153.06 + c(-0.01, 0.01),
                             rt = 243.39 + c(-60, 60), aggregationFun = "max")
#Stachydrine
eic_Stachydrine <- chromatogram(processed_Data, mz = 	144.10 + c(-0.01, 0.01),
                                rt = 60.84 + c(-60, 60), aggregationFun = "max")

#homocysteine
eic_homocysteine <- chromatogram(processed_Data, mz = 	385.13 + c(-0.01, 0.01),
                                 rt = 72.68 + c(-60, 60), aggregationFun = "max")
#Guanine
eic_Guanine <- chromatogram(processed_Data, mz = 	152.06 + c(-0.01, 0.01),
                            rt = 110.74 + c(-60, 60), aggregationFun = "max")
#Dinoxanthin
eic_Dinoxanthin <- chromatogram(processed_Data, mz = 	641.42 + c(-0.01, 0.01),
                                rt = 618.52 + c(-60, 60), aggregationFun = "max")
#Riboflavin
eic_Riboflavin <- chromatogram(processed_Data, mz = 	377.15 + c(-0.01, 0.01),
                               rt = 197.20 + c(-60, 60), aggregationFun = "max")
#Cholesterol
eic_Cholesterol <- chromatogram(processed_Data, mz = 	369.35 + c(-0.01, 0.01),
                                rt = 650.76 + c(-60, 60), aggregationFun = "max")
#Ile_Lys
eic_Ile_Lys <- chromatogram(processed_Data, mz = 	260.20 + c(-0.01, 0.01),
                            rt = 151.55 + c(-60, 60), aggregationFun = "max")

#Loliolide
eic_Loliolide <- chromatogram(processed_Data, mz = 	197.12 + c(-0.01, 0.01),
                              rt = 259.28 + c(-60, 60), aggregationFun = "max")
#myrcene
eic_Myrcene <- chromatogram(processed_Data, mz = 	159.12 + c(-0.01, 0.01),
                            rt = 388.85 + c(-60, 60), aggregationFun = "max")
#Astaxanthin
eic_Astaxanthin <- chromatogram(processed_Data, mz = 	597.39 + c(-0.01, 0.01),
                                rt = 606.32 + c(-60, 60), aggregationFun = "max")
#LUMICHROME
eic_Lumichrome <- chromatogram(processed_Data, mz = 	243.09 + c(-0.01, 0.01),
                               rt = 256.55 + c(-60, 60), aggregationFun = "max")
#Diadinoxanthin
eic_Diadinoxanthin <- chromatogram(processed_Data, mz = 	581.40 + c(-0.01, 0.01),
                                   rt = 567.38 + c(-60, 60), aggregationFun = "max")
#BETAINE
eic_Betain <- chromatogram(processed_Data, mz = 	118.09 + c(-0.01, 0.01),
                           rt = 60.21 + c(-60, 60), aggregationFun = "max")

#BiakamideC,D
eic_Biakamide <- chromatogram(processed_Data, mz = 	524.27 + c(-0.01, 0.01),
                              rt = 437.32 + c(-60, 60), aggregationFun = "max")
#Fucoxanthin
eic_Fucoxanthin <- chromatogram(processed_Data, mz = 	659.43 + c(-0.01, 0.01),
                                rt = 526.04 + c(-60, 60), aggregationFun = "max")
#Arachidonic acid 
eic_Arachidonic_acid  <- chromatogram(processed_Data, mz = 	305.25 + c(-0.01, 0.01),
                                      rt = 462.80 + c(-60, 60), aggregationFun = "max")
#Phaeophorbide
eic_Phaeophorbide  <- chromatogram(processed_Data, mz = 	593.27 + c(-0.01, 0.01),
                                   rt = 566.14 + c(-60, 60), aggregationFun = "max")
#Crassostreaxanthin
eic_Crassostreaxanthin  <- chromatogram(processed_Data, mz = 	599.41 + c(-0.01, 0.01),
                                        rt = 563.04 + c(-60, 60), aggregationFun = "max")

#Shinorine
eic_Shinorine  <- chromatogram(processed_Data, mz = 	333.13 + c(-0.01, 0.01),
                               rt = 62.44 + c(-60, 60), aggregationFun = "max")
#Pinolenic acid
eic_Pinolenic_acid  <- chromatogram(processed_Data, mz = 	161.10 + c(-0.01, 0.01),
                                    rt = 277.36 + c(-60, 60), aggregationFun = "max")
#Trimethyllysine
eic_Trimethyllysine  <- chromatogram(processed_Data, mz = 	189.16 + c(-0.01, 0.01),
                                     rt = 55.49 + c(-60, 60), aggregationFun = "max")
#Usujirene
eic_Usujirene <- chromatogram(processed_Data, mz = 	285.14 + c(-0.01, 0.01),
                              rt = 148.69 + c(-60, 60), aggregationFun = "max")
#ANTHRANILATE
eic_Anthranilate <- chromatogram(processed_Data, mz = 	120.04 + c(-0.01, 0.01),
                                 rt = 	81.32 + c(-60, 60), aggregationFun = "max")
#Succinoadenosine
eic_Succinoadenosine <- chromatogram(processed_Data, mz = 	384.12 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#Ala-Phe
eic_Ala_Phe <- chromatogram(processed_Data, mz = 	237.12 + c(-0.01, 0.01),
                            rt = 	183.02 + c(-60, 60), aggregationFun = "max")
#Diphenylphosphate
eic_Diphenylphosphate <- chromatogram(processed_Data, mz = 	251.05 + c(-0.01, 0.01),
                                      rt = 	521.71 + c(-60, 60), aggregationFun = "max")
#Octinoxate
eic_Octinoxate <- chromatogram(processed_Data, mz = 	291.19 + c(-0.01, 0.01),
                               rt = 	438.85 + c(-60, 60), aggregationFun = "max")
#Val-Trp
eic_Val_Trp <- chromatogram(processed_Data, mz = 	304.17 + c(-0.01, 0.01),
                            rt = 	209.34 + c(-60, 60), aggregationFun = "max")
#Inosine
eic_Inosine <- chromatogram(processed_Data, mz = 	269.09 + c(-0.01, 0.01),
                            rt = 	92.87 + c(-60, 60), aggregationFun = "max")
#methyl coumarin
eic_Methyl_coumarin <- chromatogram(processed_Data, mz = 	161.06 + c(-0.01, 0.01),
                                    rt = 	317.59 + c(-60, 60), aggregationFun = "max")
#Smenamide A
eic_Smenamide <- chromatogram(processed_Data, mz = 501.25 + c(-0.01, 0.01),
                              rt = 	432.08 + c(-60, 60), aggregationFun = "max")
#Ser-Ile
eic_Ser_Ile <- chromatogram(processed_Data, mz = 219.13 + c(-0.01, 0.01),
                            rt = 	153.19 + c(-60, 60), aggregationFun = "max")
#Phe-Glu
eic_Phe_Glu <- chromatogram(processed_Data, mz = 295.13 + c(-0.01, 0.01),
                            rt = 	161.04 + c(-60, 60), aggregationFun = "max")
#Petasol
eic_Petasol <- chromatogram(processed_Data, mz = 235.17 + c(-0.01, 0.01),
                            rt = 	453.49 + c(-60, 60), aggregationFun = "max")
#Butyrylcarnitine
eic_Butyrylcarnitine <- chromatogram(processed_Data, mz = 232.15 + c(-0.01, 0.01),
                                     rt = 	168.22 + c(-60, 60), aggregationFun = "max")
#palmitoyl carnitine
eic_palmitoyl_carnitine <- chromatogram(processed_Data, mz = 400.34 + c(-0.01, 0.01),
                                        rt = 	476.56 + c(-60, 60), aggregationFun = "max")
#Benzothiazolesulfonic acid
eic_Benzothiazolesulfonic_acid <- chromatogram(processed_Data, mz = 215.98 + c(-0.01, 0.01),
                                               rt = 	211.42 + c(-60, 60), aggregationFun = "max")
#LAUROYLCARNITINE
eic_Lauroylcarnitine <- chromatogram(processed_Data, mz = 344.28 + c(-0.01, 0.01),
                                     rt = 	369.72 + c(-60, 60), aggregationFun = "max")


save.image(file = file.path(env_path,
                            "DDAmyEnvironment4_gap_fillingg_EIC13_positive1.RData"))

#PDF
pdf(file.path(plot_path,
              "DDAPeak_gap_filling_positive1.pdf"),
    width = 14, height = 10)

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


####    6) export the files ####

####    6.a) export MS1 and MS2 features (GENERAL PURPOSES####

## export the individual spectra into a .mgf file (can be used for SIRIUS)
source(file.path(base_path, "customFunctions.R"))
filteredMs2Spectra <- featureSpectra(processedData, return.type = "MSpectra")
filteredMs2Spectra <- clean(filteredMs2Spectra, all = TRUE)
filteredMs2Spectra <- formatSpectraForGNPS(filteredMs2Spectra)

writeMgfData(
  filteredMs2Spectra,
  file.path(out_path,
            "DDAms2spectra_all_for_Sirius_positive1.mgf")
)

#Export peak area quantification table. 
featuresDef <- featureDefinitions(processed_Data)
featuresIntensities <- featureValues(processed_Data, value = "into")
dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
head(dataTable)
write.table(dataTable,
            file.path(out_path, "DDAfeature_table_positive1.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

####    6.b) Export MS2 features (for GNPS) ####

filteredMs2Spectra
## Select for each feature the Spectrum2 with the largest TIC.
filteredMs2Spectra_maxTic <- MSnbase::combineSpectra(
  filteredMs2Spectra,
  fcol = "feature_id",
  method = maxTic)
filteredMs2Spectra
filteredMs2Spectra_maxTic #We should observe less spectra
writeMgfData(filteredMs2Spectra_maxTic,
             file.path(out_path,
                       "DDAms2spectra_maxTic_GNPS_positive1.mgf"))

## filter data table to contain only peaks with MSMS DF[ , !(names(DF) %in% drops)]
filteredDataTable <- dataTable[which(
  dataTable$Row.names %in% filteredMs2Spectra@elementMetadata$feature_id),]
head(filteredDataTable) 
write.table(filteredDataTable,
            file.path(out_path,
                      "DDAfeature_table_xcms_onlyMS2_GNPS_positive1.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


