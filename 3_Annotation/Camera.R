##### CAMERA ####

####1) Data Loading ####

#Set the environment:
library(here)
base_path <- here()
input_path <- file.path(base_path, "1_Input_files")
env_path <- file.path(base_path, "2_Environments")
out_path <- file.path(base_path, "3_Output_files")
dir.create(env_path, recursive = TRUE, showWarnings = FALSE)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

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

#BiocManager::install("CAMERA")
library("CAMERA")

# Load the preprocessed XCMS environment generated during the data preprocessing workflow.
# The object 'mse1' is expected to be present in this environment.
load(file.path(input_path,
               "myEnvironment5_gap_filling_positive1.RData"))

#Convert your object:
xset <- as(filterMsLevel(mse1, msLevel = 1L), "xcmsSet")
xset

#Manually add the sample groups:
sampclass(xset) <- sampleData(mse1)$ESSENTIAL_sample_type
table(sampclass(xset))
head(data.frame(
  sampclass = sampclass(xset),
  metadata = sampleData(mse1)$ESSENTIAL_sample_type
))

xsa <- xsAnnotate(xset, polarity = "positive")
rm(list = setdiff(ls(), "xsa"))
gc()

#Save the progress
save.image(file = file.path(env_path,
                            "myEnvironment1_Data_Loading1.RData"))

####2) groupFWHM ####
load(file.path(env_path,
               "myEnvironment1_Data_Loading1.RData"))

#xseic_Anthranilate#xsaF <- groupFWHM(xsa, perfwhm=0.6)
xsaF <- groupFWHM(xsa, sigma = 6, perfwhm = 1) 
save.image(file = file.path(env_path,
                            "2myEnvironment_groupFWHM.RData"))

####3) findIsotopes #####
load(file.path(env_path,
               "2myEnvironment_groupFWHM.RData"))

xsaFI <- findIsotopes(
  xsaF,
  ppm = 5,
  mzabs = 0.01,
  maxcharge = 2,
  maxiso = 3,
  minfrac = 0.5,
  intval = "into"
) #Found isotopes: 22383 

#xsaFI <- findIsotopes(xsaC, maxcharge = 2, maxiso = 3, minfrac = 0.5,
#                      ppm = 5, intval = "maxo")
save.image(file = file.path(env_path,
                            "3myEnvironment_findIsotopes.RData"))

####4) groupCorr ####

xsaC <- groupCorr(
  xsaFI,
  cor_eic_th = 0.9,
  cor_exp_th = 0.8,
  pval = 0.05,
  graphMethod = "hcs",
  calcIso = TRUE,
  calcCiS = TRUE,
  calcCaS = FALSE   
)
#xsAnnotate has now 60718 groups, instead of 126 

save.image(file = file.path(env_path,
                            "4myEnvironment_groupCorr.RData"))

table(xsaC@pspectra)

####5) findAdducts ####
load(file.path(env_path,
               "4myEnvironment_groupCorr.RData")) 
rm(list = setdiff(ls(), "xsaC"))
gc()


library(BiocParallel)

param <- SnowParam(
  workers = 4,        
  type = "SOCK",      
  progressbar = TRUE
)
register(param)
#10:24
xsaFA <- findAdducts(xsaC, polarity = "positive", 
                     max_peaks = 30, multiplier = 3, ppm = 5)
cat("findAdducts finished at:", Sys.time(), "\n")

rm(list = setdiff(ls(), "xsaFA"))
gc()
save(xsaFA,
     file = file.path(env_path,
                      "xsaFA.RData"))

write.csv(getPeaklist(xsaFA),
          file = file.path(out_path,
                           "result_CAMERA.csv"))


