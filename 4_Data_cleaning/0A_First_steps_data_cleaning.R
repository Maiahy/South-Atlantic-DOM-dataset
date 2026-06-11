####Data cleaning and filtering for multivariate analysis####

#this code is based on the worked of different groups, for more information please consult:
#https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html
#https://github.com/DorresteinLaboratory/XCMS3_FeatureBasedMN/blob/master/XCMS3_Preprocessing.Rmd
#https://pietrofranceschi.github.io/Metabolomics_lectures/

#This code was written based on the recommendations provided by Broadhurst and al. 
#(DOI: https://doi.org/10.1007/s11306-018-1367-3). Please cite this article if you are using this code.
#And with the help of this website:
#https://www.davidzeleny.net/anadat-r/doku.php/en:pcoa_nmds


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("BiocParallel")

#BiocManager::install(c("xcms", "CAMERA"))

#install.packages("faahKO")
#install.packages("RColorBrewer")
#install.packages("pander")
#install.packages("magrittr")
#install.packages("pheatmap")
#install.packages("SummarizedExperiment")


#Loading of the packages:
library(xcms)
library(RColorBrewer)
library(pander)
library(magrittr)
library(SummarizedExperiment)
library(ggplot2)
library(tidyverse)
library(readxl)
library(tidyr)
library(dplyr)
library(BiocParallel)


register(SnowParam(7))

#Environment:
library(here)
base_path <- here()
input_path <- file.path(base_path, "1_Input_files")
env_path <- file.path(base_path, "2_Environments")
plot_path <- file.path(base_path, "3_Quality_plots")
out_path <- file.path(base_path, "4_Output_files")
dir.create(env_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)


##### Data import #### 


# feature table
feature_table <- read.table(
  file.path(input_path,
            "feature_table_MS1_positive1.txt"),
  header = TRUE
)
colnames(feature_table)
colnames(feature_table)[which(names(feature_table) == "Row.names")] <- "featureId"

feature_table <- feature_table [, ! colnames(feature_table) %in% 
                                  c("X","File_name")]

colnames(feature_table)

#read the metadata of the samples
library(janitor)
metadata_samples <- read.table(
  file.path(input_path,
            "Metadata_final_with_operators_withoutAcquireX.txt"),
  sep = "\t"
)
colnames(metadata_samples)
metadata_samples <- metadata_samples %>%
  mutate(ESSENTIAL_sample_type = case_when(
    ESSENTIAL_sample_id == "InterbatchQC" ~ "REF",
    str_detect(ESSENTIAL_sample_id, "AcquireX") ~ "AcquireX",
    TRUE ~ ESSENTIAL_sample_type
  ))

metadata_samples %>%
  filter(str_detect(ESSENTIAL_sample_type, "QC")) %>%   # ou == "QC" si exact
  distinct(ESSENTIAL_sample_id)
unique(metadata_samples$ESSENTIAL_sample_type)

library(ggplot2)
library(tidyverse)
library(readxl)
library(tidyr)
library(dplyr)

DM <- feature_table %>% 
  column_to_rownames("featureId") %>% 
  t(.) %>% 
  as_tibble(rownames = "file_name")

#Remove the X at the begenning of some lines:
DM$file_name = gsub("^X", "", DM$file_name)

#Remove all the lines that are not startting by a date:
DM = DM[grepl("^\\d{8}", DM$file_name),]

#to save the environment (useful if very long process):
save.image(file=  file.path(env_path,'myEnvironment1_data_import1.RData'))
rm(feature_table)

colnames(metadata_samples)


#####1) Normalization based on the volume####


metadata_samples$Volume_normalization_factor[is.na(metadata_samples$Volume_normalization_factor)] <- 1
metadata_samples <- metadata_samples[, c(setdiff(names(metadata_samples), "Volume_normalization_factor"), "Volume_normalization_factor")]


# Select specific columns from the dataframe
metadata_samples2 <- metadata_samples %>%
  select(Volume_normalization_factor, file_name)


#Left_joint the both tables in order to obtain a working table:
Working_table = left_join(DM, metadata_samples2, by = "file_name")

#A normalization based n the volume can be performed in 2 steps:
#1) Calculatation of a dilution factore with DF = V sample/ Volume reference. Here we consider the reference volume to be 2L
#2) Normalization of the intensities based on the dilution factor: divide the intensity of the metabolite by the dilution factor.
#Here we have already calculated DF in the excel file we then just have to divide the entire table by the values of DF contained in the metadata.
Working_table <- as.data.frame(Working_table)
# Number of columns in the dataframe
num_cols <- ncol(Working_table)
num_cols <- as.numeric(num_cols)

# Loop through each row
for (i in 1:nrow(Working_table)) {
  # Get the value of the last column in the current row
  divisor <- as.numeric(Working_table[i, num_cols])
  
  # Convert other columns (except the first and last) to numeric
  Working_table[i, 2:(num_cols-1)] <- as.numeric(Working_table[i, 2:(num_cols-1)])
  
  # Divide all values in the row (except the first and last columns) by the divisor
  Working_table[i, 2:(num_cols-1)] <- Working_table[i, 2:(num_cols-1)] / divisor
}


# Output the modified dataframe
print(Working_table)

#to save the environment (useful if very long process):
save.image(file= file.path(env_path,'myEnvironment1_normalisation1.RData'))

#Export comparison tables to see if it was working
colnames(DM)
DM_before_normalization = DM %>%
  select(file_name, FT00036, FT00829)
write.table(
  DM_before_normalization,
  file = file.path(out_path,
                   "DM_before_normalization.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


unique(Working_table$VOLUME_normalization_factor)
unique(metadata_samples$VOLUME_normalization_factor)

Working_table$VOLUME_normalization_factor
Working_table_after_normalization = Working_table %>%
  select(file_name, FT00036, FT00829, Volume_normalization_factor)
write.table(
  Working_table_after_normalization,
  file = file.path(
    out_path,
    "Working_table_after_normalization.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
) # Compare the manually calculated values with the values generated in R:


Working_table <- Working_table[, -ncol(Working_table)] #remove the last column with normalizaion factor
write.table(
  Working_table,
  file = file.path(
    out_path,
    "feature_table_normalised_volume.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
save.image(file=  file.path(env_path,'myEnvironment2_volume_normalization.RData'))


#####2) Blank filtering ####

#This method is inspired from McParland and al. 
#"Seasonal patterns of DOM molecules are linked to microbial functions in the oligotrophic ocean"
#https://doi.org/10.1128/msystems.01540-25
#Please consult their article

 
#1) calculation of the average value in blank and in samples for each 
#feature
DM2 <- Working_table%>% 
  left_join(metadata_samples, by = "file_name")
colnames(DM2)

DM2 %>%
  filter(ESSENTIAL_sample_type == "Blank") %>%
  pull(ESSENTIAL_sample_id) %>%
  unique()

DM2 %>%
  filter(ESSENTIAL_sample_type == "Blank") %>%
  pull(ESSENTIAL_sample_id) %>%
  table()


DM2 <- DM2 %>%
  mutate(blank_category = case_when(
    ESSENTIAL_sample_type == "Blank" & grepl("^SAMEA", ESSENTIAL_sample_id) ~ "onboard_blank",
    ESSENTIAL_sample_type == "Blank" & grepl("ProcessBlank", ESSENTIAL_sample_id) ~ "process_blank",
    ESSENTIAL_sample_type == "Blank" & grepl("MeOHBlank", ESSENTIAL_sample_id) ~ "meoh_blank",
    TRUE ~ NA_character_
  ))

DM2 %>%
  filter(ESSENTIAL_sample_type == "Blank") %>%
  pull(blank_category) %>%
  table()

DM2_log <- DM2 %>%
  mutate(across(starts_with("FT"), ~ log2(.x + 1)))

samples_average = DM2_log %>%
  filter(ESSENTIAL_sample_type == "Sample") %>%
  summarise(across(starts_with("FT"),
                   mean, na.rm = TRUE)) %>%
  mutate(across(starts_with("FT"),
                ~ ifelse(is.nan(.x), 0, .x)))

blanks_average_by_type <- DM2_log %>%
  filter(ESSENTIAL_sample_type == "Blank") %>%
  group_by(blank_category) %>%
  summarise(across(starts_with("FT"),
                   mean, na.rm = TRUE)) %>%
  mutate(across(starts_with("FT"),
                ~ ifelse(is.nan(.x), 0, .x)))


blanks_max <- blanks_average_by_type %>%
  summarise(across(starts_with("FT"),
                   ~ max(.x, na.rm = TRUE)))


blanks_max[is.infinite(as.matrix(blanks_max))] <- 0


all(colnames(samples_average) == colnames(blanks_max))
df_md <- tibble(
  feature = colnames(samples_average),
  sample = as.numeric(samples_average[1,]),
  blank = as.numeric(blanks_max[1,])
) %>%
  mutate(
    diff = sample - blank,
    mean = (sample + blank) / 2
  )

library(ggplot2)

ggplot(df_md, aes(x = mean, y = diff)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, color = "blue") +
  labs(
    x = "Mean intensity",
    y = "Difference (Sample - Blank)",
    title = "Mean-Difference Plot"
  ) +
  theme_minimal()

artifact <- df_md %>%
  filter(mean < 15 & diff > 10)

artifact_check <- tibble(
  feature = artifact$feature,
  sample = artifact$sample,
  blank = artifact$blank
)

head(artifact_check)
summary(artifact_check$blank)
artifact_check %>%
  summarise(
    mean_sample = mean(sample),
    median_sample = median(sample)
  )


quantiles <- quantile(df_md$mean,
                      probs = c(0.2, 0.4, 0.6, 0.8, 1),
                      na.rm = TRUE)

quantiles
df_md <- df_md %>%
  mutate(bin = cut(mean,
                   breaks = c(-Inf, quantiles),
                   labels = c("bin1","bin2","bin3","bin4","bin5")))

ggplot(df_md, aes(mean, diff, color = bin)) +
  geom_point(alpha = 0.3) +
  theme_minimal()


thresholds <- df_md %>%
  filter(diff < 0) %>%
  group_by(bin) %>%
  summarise(threshold = abs(quantile(diff, 0.25, na.rm = TRUE)))

df_md <- df_md %>%
  left_join(thresholds, by = "bin")

df_md <- df_md %>%
  mutate(keep = diff > threshold)


ggplot(df_md, aes(mean, diff, color = keep)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0) +
  theme_minimal()


df_filtered <- df_md %>%
  filter(keep == TRUE)
features_to_keep <- df_filtered$feature

DM3 <- DM2 %>%
  select(file_name, starts_with("FT")) #99207 features

DM3_filtered <- DM3 %>%
  select(file_name, all_of(features_to_keep)) # 23585 features

save.image(file=  file.path(env_path,'myEnvironment3_blanksubstraction.RData'))
load(file.path(env_path,'myEnvironment3_blanksubstraction.RData'))

######b) Figure####
quantiles <- quantile(df_md$mean,
                      probs = c(0.2, 0.4, 0.6, 0.8, 1),
                      na.rm = TRUE)

df_md <- df_md %>%
  mutate(bin = cut(mean,
                   breaks = c(-Inf, quantiles),
                   labels = c("0–20%", "20–40%", "40–60%", "60–80%", "80–100%"),
                   include.lowest = TRUE))
p_bins <- ggplot(df_md, aes(x = mean, y = diff, color = bin)) +
  geom_point(
    alpha = 0.25,   # un peu plus léger = mieux visuellement
    size = 1
  ) +
  
  scale_color_manual(
    values = c("#c7c72e", "#61b347", "#4bbcbe", "#579ad1", "#b97cb4"),
    drop = FALSE   # évite de perdre une classe si vide
  ) +
  
  labs(
    x = "Mean intensity (log2)",
    y = "Difference (Sample - Blank)",
    color = "Intensity percentiles"
  ) +
  
  theme_classic(base_size = 14) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "plain"),
    legend.position = "right",
    legend.title = element_text(face = "plain"),
    legend.key.size = unit(0.4, "cm"),
    aspect.ratio = 1
  )

p_keep <- ggplot(df_md, aes(x = mean, y = diff, color = keep)) +
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c("TRUE" = "#579ad1", "FALSE" = "#df6765"),
    labels = c("FALSE" = "Removed", "TRUE" = "Retained")
  ) +
  
  labs(
    x = "Mean intensity (log2)",
    y = "Difference (Sample - Blank)",
    color = "Filtering"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "plain"),
    legend.position = "right",
    aspect.ratio = 1
  )

library(patchwork)
final_plot <- p_bins + p_keep +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(plot_path,"Blank_filtering_panel.pdf"),
       plot = final_plot,
       width = 12,
       height = 6)

######c) Check the results####


df_md %>%
  group_by(bin, keep) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin) %>%
  mutate(prop = n / sum(n))


thresholds

#####3)	Replicate filtration  ####

rm(list = setdiff(ls(), c("DM3_filtered", "metadata_samples")))

DM4 <- DM3_filtered %>% 
  left_join(metadata_samples, by = "file_name")
tail(colnames(DM4), 100)

unique(DM4$ESSENTIAL_sample_type)

DM4 <- DM4 %>%
  filter(ESSENTIAL_sample_type == "Sample")  

library(tidyr)

DM4_long <- DM4 %>%
  pivot_longer(
    cols = starts_with("FT"),
    names_to = "feature",
    values_to = "intensity"
  )

DM4_long <- DM4_long %>%
  mutate(detected = intensity > 0,
         detected = replace_na(detected, FALSE))


replicate_summary <- DM4_long %>%
  group_by(ESSENTIAL_sample_id, feature) %>%
  summarise(
    n_detected = sum(detected),
    .groups = "drop"
  )

replicate_summary <- replicate_summary %>%
  filter(n_detected > 0)

feature_stability <- replicate_summary %>%
  mutate(valid_group = n_detected >= 2) %>%
  group_by(feature) %>%
  summarise(
    n_groups = n(),
    n_valid = sum(valid_group),
    proportion_valid = n_valid / n_groups,
    .groups = "drop"
  )

sum(feature_stability$proportion_valid >= 0.7)  # gardées 10403 features
sum(feature_stability$proportion_valid < 0.7)   # supprimées 15817

features_to_keep <- feature_stability %>%
  filter(proportion_valid >= 0.7) %>%
  pull(feature)


DM4_filtered <- DM3_filtered %>%
  select(file_name, all_of(features_to_keep)) # 10127 features
save.image(file =  file.path(env_path,'myEnvironment4_Replicatefiltration.RData'))

#####4)	Detection rate ####
rm(list = setdiff(ls(), c("DM4_filtered", "metadata_samples")))

DM5 <- DM4_filtered %>% 
  left_join(metadata_samples, by = "file_name")

QC_only <- DM5 %>%
  filter(ESSENTIAL_sample_type == "REF")

ft_columns <- grep("^FT", names(QC_only), value = TRUE)

QC_only2 <- QC_only[, ft_columns]

detection_rate <- colMeans(!is.na(QC_only2))

Detection_df <- data.frame(
  feature = names(detection_rate),
  detection_rate = detection_rate
)

not_detected <- Detection_df %>%
  filter(detection_rate == 0) %>%
  pull(feature)
length(not_detected)

low_detection_features <- Detection_df %>%
  filter(detection_rate < 0.7) %>%
  pull(feature)


Detection_df %>%
  summarise(
    total = n(),
    low_detect = sum(detection_rate < 0.7),
    high_detect = sum(detection_rate >= 0.7),
    prop_low = low_detect / total
  )

ggplot(Detection_df, aes(x = detection_rate)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = 0.7, color = "red", linewidth = 0.5  ) +
  theme_minimal()

p <- ggplot(Detection_df, aes(x = detection_rate)) +
  geom_density(
    fill = "#579ad1",
    color = "#579ad1",
    alpha = 0.4,
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0.7,
    color = "black",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "Detection frequency in pooled QC samples",
    y = "Feature density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "plain")
  )

ggsave(file.path(plot_path,"Detection_rate.pdf"), plot = p,
       width = 7, height = 4)

high_detection_features <- Detection_df %>%
  filter(detection_rate >= 0.7) %>%
  pull(feature)

DM5_filtered <- DM4_filtered %>%
  select(file_name, all_of(high_detection_features))

save.image(file=  file.path(env_path,'myEnvironment5_DetectionRate.RData'))


######b) an other figure ####
zero_QC_features <- Detection_df %>%
  filter(detection_rate == 0) %>%
  pull(feature)


DM_samples <- DM5 %>%
  filter(ESSENTIAL_sample_type == "Sample") %>%
  select(file_name, all_of(zero_QC_features))

library(tidyr)

DM_samples_long <- DM_samples %>%
  pivot_longer(
    cols = -file_name,
    names_to = "feature",
    values_to = "intensity"
  )

library(ggplot2)

p_zeroQC <- ggplot(DM_samples_long, aes(x = intensity)) +
  geom_density(
    fill = "#df6765",
    color = "#df6765",
    alpha = 0.4,
    linewidth = 0.2
  ) +
  labs(
    x = "Signal intensity (a.u.)",
    y = "Density",
    title = "Features not detected in QC samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )

ggsave(file.path(plot_path,"Zero_QC_features_density.pdf"),
       plot = p_zeroQC,
       width = 7,
       height = 4)

#####5) Precision ####
rm(list = setdiff(ls(), c("DM5_filtered", "metadata_samples", "feature_annotation")))

#ISO 5725 define the precision as the “ measure of dispersion of the distribution 
#of [independent] test results”, where, “independent test results are obtained with
#the same method on identical test items in the same laboratory by the same operator 
#using the same equipment within short intervals of time”.

#This could therefore be measured with our pooled QCs, that are injected multiple times
#in our system.

#The convention dictate that the RSD of one feature across the different QC should 
#not be higher than 20% for a LC-MS system and 30% for a GC-MS system.

#Columns of metadata to keep:
colnames(metadata_samples)

metadata_samples2 = metadata_samples %>% 
  select(file_name, ESSENTIAL_sample_id, biological_replicate)
colnames(metadata_samples2)

#Fusion the information about the sample + the biological replicate to see the unique sample measured
metadata_samples2$sample <- paste(metadata_samples2$ESSENTIAL_sample_id, metadata_samples2$biological_replicate)
metadata_samples2 = metadata_samples2 %>% 
  select(file_name, sample)
library(stringr)
metadata_samples2 <- as.data.frame(
  apply(metadata_samples2,2, function(x) gsub("\\s+", "", x)))

DM2 <- DM5_filtered %>% 
  left_join(metadata_samples2, by = "file_name")
colnames(DM2)

DM2$file_name <- NULL

DM2 <- DM2 %>%
  mutate(
    sample = dplyr::case_when(
      sample %in% c("InterbatchQCR01", "InterbatchQCR02") ~ "InterbatchQC",
      TRUE ~ sample
    )
  )

#Calculate median of intensisty for each metabolite into each sample accross the replicates:
feature_median_per_sample = DM2 %>%
  group_by(sample) %>%
  summarise(across(starts_with("FT"), median, na.rm = TRUE))

#Calculate MAD per sample:
feature_mad_per_sample = DM2 %>%
  group_by(sample) %>%
  summarise(across(starts_with("FT"), mad, na.rm = TRUE))


#Remove all the lines that are not concerning the QC in the both tables:
feature_median_QC = subset(feature_median_per_sample, sample == "InterbatchQC")
feature_mad_QC = subset(feature_mad_per_sample, sample == "InterbatchQC")

colnames(feature_mad_QC)
class(feature_mad_QC$FT00002)
class(feature_median_QC$FT00002)

#RSD_QC = ((feature_mad_QC * 1.4826)/ feature_median_QC) *100%
RSD_QC <- cbind(
  feature_mad_QC[1],
  round(
    (1.4826 * feature_mad_QC[-1]) / 
      ifelse(feature_median_QC[-1] == 0, NA, feature_median_QC[-1]) * 100,
    1
  )
)

###
###
#PLOT#

# Extraire les noms des features (première colonne de RSD_QC)
feature_names <- RSD_QC[, 1]

# Extraire les RSDs (colonnes 2 à ncol(RSD_QC))
RSD_values <- as.numeric(RSD_QC[, -1])

# Créer un data frame pour ggplot
RSD_df <- data.frame(
  Feature = feature_names,
  RSD = RSD_values)


library(ggplot2)
RSD_df_clean <- RSD_df %>%
  filter(!is.na(RSD))

pdf(
  file.path(
    plot_path,
    "RSD_before_SERRF_and_filtering_50.pdf"
  ),
  width = 14,
  height = 10
)
ggplot(RSD_df_clean, aes(x = RSD)) +
  
  # histogram
  geom_histogram(
    bins = 50,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  labs(
    x = "RSD in QC samples (%)",
    y = "Number of features"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "plain"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

dev.off()

RSD_summary <- table(cut(RSD_df$RSD,
                         breaks = c(-Inf, 20, 30, 50, Inf),
                         labels = c("0-20%", "20-30%", "30-50%", ">50%")))
print(RSD_summary)



#0-20% 20-30% 30-50%   >50% 
 # 2    549   2697   4642 

sum(is.na(RSD_values)) #0

##
##

#2) Let#s now collect the columns for which the value is superior to 50
#we then consider that the variation in the QC is to high and that the features cannot being studied,
#even after batch correction.
selected_columns <- colnames(RSD_QC)[-1][
  is.na(RSD_values) | RSD_values > 50
]

head(selected_columns)
length(selected_columns)

selected_columns = as.character(selected_columns)
selected_columns = as.list(selected_columns)

#3) Let`s remove all the columns corresponding to some noise in the data table
DM6_filtered= DM5_filtered[ , -which(names(DM5_filtered) %in% selected_columns)] #3248 features


#Save the data table for further analysis:
write.table(
  DM6_filtered,
  file = file.path(
    out_path,
    "0_final_feature_table_univariate_analysis_before_batch_correction_RSD50.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#Save your work
save.image(file=  file.path(env_path,'myEnvironment6_RSD1.RData'))  

unique(metadata_samples$ESSENTIAL_sample_id)


######b)RSD of the other QCS before SERRF correction #####


batch_qc_ids <- c("Batch1QC","Batch2QC","Batch3QC","Batch4QC")

# Ajouter metadata
DM_QC <- DM6_filtered %>%
  left_join(metadata_samples, by = "file_name") %>%
  filter(ESSENTIAL_sample_id %in% batch_qc_ids)

# Fonction pour calculer RSD
compute_rsd <- function(df) {
  
  median_vals <- df %>%
    summarise(across(starts_with("FT"), median, na.rm = TRUE))
  
  mad_vals <- df %>%
    summarise(across(starts_with("FT"), mad, na.rm = TRUE))
  
  rsd <- (1.4826 * mad_vals) /
    ifelse(median_vals == 0, NA, median_vals) * 100
  
  return(as.numeric(rsd))
}

# Appliquer par batch QC
RSD_list <- DM_QC %>%
  group_split(ESSENTIAL_sample_id)

names(RSD_list) <- unique(DM_QC$ESSENTIAL_sample_id)

RSD_results <- lapply(RSD_list, compute_rsd)

RSD_df <- do.call(cbind, RSD_results)

colnames(RSD_df) <- names(RSD_results)

RSD_df <- as.data.frame(RSD_df)

RSD_summary_list <- lapply(RSD_df, function(x) {
  table(cut(x,
            breaks = c(-Inf, 20, 30, 50, Inf),
            labels = c("0-20%", "20-30%", "30-50%", ">50%")))
})
RSD_summary_list



#####6) Missing values imputation ####
load(file.path(env_path,'myEnvironment6_RSD1.RData'))
rm(list = setdiff(ls(), c("DM6_filtered", "metadata_samples", "feature_annotation")))


DM6_long = DM6_filtered %>%
  pivot_longer(
    cols = starts_with("FT"),
    names_to = "feature",
    values_to = "intensity"
  ) 

#Save the exact position of the NAs
na_long <- DM6_filtered %>%
  pivot_longer(-file_name, names_to = "feature", values_to = "value") %>%
  filter(is.na(value)) %>%
  select(file_name, feature)
write.table(
  na_long,
  file = file.path(
    out_path,
    "positions_na.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

myimputer <- function(v) {
  if (sum(is.na(v)) == 0) {  
    # If there are no NAs, leave the vector unchanged
    return(v)
  } else {
    napos <- which(is.na(v))  # Find the positions of the NAs in the vector
    min_val <- min(v, na.rm = TRUE)  # Calculate the minimum non-NA value
    
    # Determine the range for random values based on the minimum value
    if (min_val > 2E5) {
      newval <- runif(length(napos), 0, 2E5)  # Random values between 0 and 1E4
    } else {
      newval <- runif(length(napos), 0, min_val)  # Random values between 0 and min_val
    }
    
    out <- v
    out[napos] <- newval  # Replace NAs with the new random values
    return(out)
  }
}

DM7_imputed <- DM6_filtered %>% 
  mutate(across(starts_with("FT"), ~ myimputer(.x)))


#Visualisation of the result:
DM7_long <- DM7_imputed %>%
  pivot_longer(
    cols = starts_with("FT"),
    names_to = "feature",
    values_to = "intensity"
  )

imputed_values <- left_join(na_long, DM7_long, by = c("file_name", "feature"))  %>%
  mutate(was_NA = "true")


observed_values <- DM7_long %>%
  anti_join(na_long, by = c("file_name", "feature")) %>%
  mutate(was_NA = "false")

final_df <- bind_rows(imputed_values, observed_values)
any(is.na(final_df))

ggplot(final_df, aes(x = intensity, fill = was_NA, color = was_NA)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("false" = "#579ad1", "true" = "#df6765")) +
  scale_color_manual(values = c("false" = "#579ad1", "true" = "#df6765")) +
  #scale_x_log10() +
  xlim(0, 4E5) +
  labs(
    x = "Intensity",
    y = "Density"
  ) +
  theme_classic(base_size = 14)


ggplot(final_df, aes(x = intensity, fill = was_NA, color = was_NA)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("false" = "#579ad1", "true" = "#df6765")) +
  scale_color_manual(values = c("false" = "#579ad1", "true" = "#df6765")) +
  #scale_x_log10() +
  xlim(0, 5E5) +
  labs(
    x = "Intensity",
    y = "Density"
  ) +
  theme_classic(base_size = 14)


p = ggplot(final_df, aes(x = intensity, fill = was_NA, color = was_NA)) +
  geom_density(alpha = 0.3) +
  
  scale_fill_manual(
    values = c("false" = "#579ad1", "true" = "#df6765"),
    labels = c("false" = "Observed values", 
               "true"  = "Imputed values")
  ) +
  
  scale_color_manual(
    values = c("false" = "#579ad1", "true" = "#df6765"),
    labels = c("false" = "Observed values", 
               "true"  = "Imputed values")
  ) +
  
  xlim(0, 4E4) +
  
  labs(
    x = "Signal intensity (a.u.)",
    y = "Density",
    fill = "Data type",
    color = "Data type"
  ) +
  
  theme_classic(base_size = 14)

ggsave(file.path(plot_path,"Imputation_method2.pdf"), plot = p,
       width = 7, height = 4)


#Save the data table for further analysis:
write.table(
  DM7_imputed,
  file = file.path(
    out_path,
    "0_feature_table_after_imputation_2E5.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

#Save your work
save.image(file=  file.path(env_path,'myEnvironment7_imputation2E5.RData'))

#####7) SERRF export ####

###Attention! It will most probably be running for several hours...
load(file.path(env_path,'myEnvironment7_imputation2E5.RData'))
rm(list = setdiff(ls(), c("DM7_imputed", "metadata_samples", "feature_annotation")))


# 1. Filtrer les types utiles

colnames(metadata_samples)
unique(metadata_samples$batch_identifier)
meta <- metadata_samples %>%
  filter(
    ESSENTIAL_sample_type %in% c("Sample", "REF", "QC"),   # types utiles
    !batch_identifier %in% c("B1", "B2", "B3", "B4", "B5")  # exclusion DDA !
  ) %>%
  select(
    file_name, batch_identifier,
    ESSENTIAL_sample_type, ESSENTIAL_sample_id, Injection_order
  )


# 2. Mapping vers les types SERRF (aucun validate)
meta$sampleType <- NA_character_

meta$sampleType[meta$ESSENTIAL_sample_type == "Sample"] <- "sample"   # échantillons réels
meta$sampleType[meta$ESSENTIAL_sample_type == "REF"]    <- "qc"       # REF = QC de calibration internes
meta$sampleType[meta$ESSENTIAL_sample_type == "QC"]     <- "sample"   # QC globaux = traités comme échantillons

# 3. Ordonner selon l’ordre d’injection
meta$Injection_order <- as.numeric(as.character(meta$Injection_order))
meta <- meta[order(meta$Injection_order), ]

# 4. Restreindre DM_i4 aux fichiers présents dans meta
meta <- meta[meta$file_name %in% DM7_imputed$file_name, ]
ft_cols <- grep("^FT", names(DM7_imputed), value = TRUE)
DM_keep <- DM7_imputed[match(meta$file_name, DM7_imputed$file_name), c("file_name", ft_cols)]

# 5. Transposer : features en lignes, samples en colonnes
mat <- t(as.matrix(DM_keep[ , -1, drop = FALSE]))
feat_ids <- seq_len(nrow(mat))
feat_labels <- rownames(mat)

# 6. Construire les 4 lignes d’en-tête attendues par SERRF
hdr1 <- c("", "batch",      as.character(meta$batch_identifier))
hdr2 <- c("", "sampleType", as.character(meta$sampleType))
hdr3 <- c("", "time",       as.character(meta$Injection_order))
hdr4 <- c("No", "label",    as.character(meta$file_name))

# 7. Combiner tout dans la table finale
serrf <- rbind(hdr1, hdr2, hdr3, hdr4,
               cbind(as.character(feat_ids), feat_labels, apply(mat, 2, as.character)))
colnames(serrf) <- NULL

# 8. Écrire le fichier CSV final
write.table(
  serrf,
  file = file.path(
    out_path,
    "for_SERRF.csv"
  ),
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

#####8) SERRF processing ####

#To run SERRF in local R:
#Open the code "app" (= Open the app.R file in RStudio.)
#Then in RStudio, click Run App button.
#Click Open in Browser. 

#At the end of the process you should get several files in your working directory,
#including the normalized data and some graphs showing the quality of your data.



