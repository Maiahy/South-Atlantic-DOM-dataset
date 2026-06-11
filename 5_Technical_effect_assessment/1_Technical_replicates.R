#PCA in R


####1) Preparation of the data ####

#Environment:
library(here)
base_path <- here()

preprocessing_path <- file.path(
  base_path,
  "2_Data_preprocessing",
  "2_Output_files"
)
# Expected objects:
# DM_i5_filtered
# metadata_samples
# feature_annotation

input_path  <- file.path(base_path, "6_Data_cleaning", "3_Files_after_cleaning")
output_path <- file.path(base_path, "7_Data_validation")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
#Loading of the packages:

library(pander)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(tidyr)
library(dplyr)
library("variancePartition")
library(janitor)

load(file.path(
  input_path,
  "myEnvironment2_RSD_after_SERRF.RData"
))


Metadata_variance_partition <- metadata_samples %>%
  filter(ESSENTIAL_sample_type == "Sample")

#Depth should be a fixed effect = continuous variable
convert_depth <- function(x) {
  if (grepl("-", x)) {
    bounds <- strsplit(x, "-")[[1]]
    return(mean(as.numeric(bounds)))
  } else {
    return(as.numeric(x))
  }
}

Metadata_variance_partition$depth_numeric <- 
  sapply(Metadata_variance_partition$ESSENTIAL_elevation_nominal_median_m, convert_depth)
Metadata_variance_partition$depth_numeric <- as.numeric(Metadata_variance_partition$depth_numeric)
Metadata_variance_partition$ESSENTIAL_event_latitude_N_midpoint <- as.numeric(Metadata_variance_partition$ESSENTIAL_event_latitude_N_midpoint)
Metadata_variance_partition$ESSENTIAL_event_longitude_E_midpoint <- as.numeric(Metadata_variance_partition$ESSENTIAL_event_longitude_E_midpoint)
Metadata_variance_partition$CTD_salinity_psu_median <- as.numeric(Metadata_variance_partition$CTD_salinity_psu_median)
Metadata_variance_partition$CTD_temperature_degC_median <- as.numeric(Metadata_variance_partition$CTD_temperature_degC_median)
Metadata_variance_partition$CTD_fluorescence_mgm3_median <- as.numeric(Metadata_variance_partition$CTD_fluorescence_mgm3_median)
Metadata_variance_partition$CTD_oxygen_percent_median <- as.numeric(Metadata_variance_partition$CTD_oxygen_percent_median)
#Create a factor for evaporation status
Metadata_variance_partition$Evaporation_status <- factor(
  Metadata_variance_partition$Evaporation_status,
  levels = c("minimal", "low", "moderate", "high"),
  ordered = TRUE
)

#Log transform the data
sample_names <- DM_i5_filtered[, 1, drop = FALSE]
intensity_data <- DM_i5_filtered[, -1]
log_transformed_data <- log2(intensity_data + 1)
DM_i5 <- cbind(sample_names, log_transformed_data)

# Extract sample names
library(mdatools)
sample_names <- DM_i5[, 1, drop = FALSE]
intensity_data <- DM_i5[, -1]
intensity_data <- as.matrix(intensity_data)
feature_table = cbind(sample_names, intensity_data)

####2) PCoA #####
metadata = Metadata_variance_partition 

feature_table <- t(feature_table)
colnames(feature_table) <- as.character(feature_table[1, ])
feature_table <- feature_table[-1, ]
feature_table <- as.data.frame(feature_table)
feature_table <- tibble::rownames_to_column(feature_table, "FeatureId")

#Load RT and filtr features
RT_infos <- read.table(
  file.path(
    preprocessing_path,
    "feature_table_MS1_positive1.txt"
  ),
  sep = "\t",
  header = TRUE,
  row.names = 1
) %>%
  tibble::rownames_to_column("FeatureId")
selected_features <- RT_infos$FeatureId[RT_infos$rtmed < 120 | RT_infos$rtmed > 500]
feature_table_filtered <- feature_table[!feature_table$FeatureId %in% selected_features, ]
feature_table_filtered <- t(feature_table_filtered)
feature_table_filtered <- row_to_names(as.data.frame(feature_table_filtered), 1)
feature_table_filtered <- tibble::rownames_to_column(feature_table_filtered, "file_name")

#Merge and filter
pca_input <- feature_table_filtered %>%
  left_join(metadata, by = "file_name") %>%
  filter(ESSENTIAL_sample_type == "Sample")



colnames(metadata)
depth_table = as.data.frame(table(metadata$depth_numeric))
pca_input_surface <- pca_input 
#pca_input_surface <- pca_input %>%
 # filter(!is.na(Evaporation_status),
     #    !is.na(CTD_temperature_degC_median),
      #   !is.na(depth_numeric))

#Matrix of features only
pca_matrix <- pca_input_surface %>%
  select(starts_with("FT")) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

#install.packages("vegan")
library(vegan)

bray_dist <- vegdist(pca_matrix, method = "bray")

pcoa_res <- cmdscale(bray_dist, eig = TRUE, add = TRUE)
#pcoa_res <- wcmdscale(bray_dist, eig = TRUE, add = TRUE)

barplot(pcoa_res$eig)
eig_vals <- pcoa_res$eig
eig_vals <- eig_vals[eig_vals > 0]
variance_explained <- eig_vals / sum(eig_vals)
pc1_var <- round(variance_explained[1] * 100, 1)
pc2_var <- round(variance_explained[2] * 100, 1)

pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")


pcoa_df$sample_id <- pca_input_surface$ESSENTIAL_sample_id
pcoa_df$file_name <- pca_input_surface$file_name
pcoa_df$station <- pca_input_surface$ESSENTIAL_station
pcoa_df$topical_study <- pca_input_surface$project_topical_study_label
library(plotly)
pcoa_df$technical_replicate <- pca_input_surface$technical_replicate


pcoa_df$technical_replicate <- factor(
  pcoa_df$technical_replicate,
  levels = c("R1", "R2", "R3")
)

pcoa_df <- pcoa_df %>%
  arrange(sample_id, technical_replicate)
library(plotly)
colnames(metadata)

######segemnts to centroid#####

library(dplyr)

centroids <- pcoa_df %>%
  group_by(sample_id) %>%
  summarise(
    cx = mean(PCoA1),
    cy = mean(PCoA2),
    .groups = "drop"
  )

pcoa_df2 <- left_join(pcoa_df, centroids, by = "sample_id")

library(plotly)
library(viridis)





##### Is intra sample > to inter sample ####

# Distance intra-échantillon
intra_dist <- pcoa_df2 %>%
  group_by(sample_id) %>%
  summarise(
    intra = mean(sqrt((PCoA1 - cx)^2 + (PCoA2 - cy)^2))
  )

# Distance inter-échantillon (centroïdes entre eux)
centroids_dist <- dist(centroids[, c("cx", "cy")])
inter_dist <- mean(centroids_dist)

mean(intra_dist$intra)
inter_dist

mean_intra <- mean(intra_dist$intra)
ratio <- mean_intra / inter_dist

mean_intra
inter_dist
ratio

#adonis2(bray_dist ~ ESSENTIAL_sample_id, data = pca_input_surface)

#Where are my outliers? 
intra_dist %>%
  ggplot(aes(x = intra)) +
  geom_histogram(bins = 30)

#Boxplot rule to remove the outliers:
Q3 <- quantile(intra_dist$intra, 0.75)
IQR_val <- IQR(intra_dist$intra)

threshold <- Q3 + 1.5 * IQR_val
threshold

bad_samples <- intra_dist %>%
  filter(intra > threshold)

bad_samples #34 samples are possible outliers... 8.7% d’outliers

bad_samples_info <- bad_samples %>%
  left_join(
    pca_input_surface,
    by = c("sample_id" = "ESSENTIAL_sample_id")
  )


library(writexl)

write_xlsx(
  bad_samples,
  path = file.path(
    output_path,
    "Outliers.xlsx"
  )
)


colnames(metadata)
bad_samples_stations = bad_samples_info %>%
  select(sample_id, intra, ESSENTIAL_station_label, project_topical_study_label, Evaporation_status, Code_inlab_operator )

table(bad_samples_stations$Evaporation_status)
table(bad_samples_stations$Code_inlab_operator)
unique(bad_samples_stations$ESSENTIAL_station_label)
table(bad_samples_stations$project_topical_study_label)

n_outliers_per_evoparion = pca_input_surface %>%
  mutate(is_outlier = ESSENTIAL_sample_id %in% bad_samples$sample_id) %>%
  group_by(Evaporation_status) %>%
  summarise(
    n_total = n(),
    n_outliers = sum(is_outlier),
    ratio = n_outliers / n_total
  )

n_outliers_per_operator = pca_input_surface %>%
  mutate(is_outlier = ESSENTIAL_sample_id %in% bad_samples$sample_id) %>%
  group_by(Code_inlab_operator) %>%
  summarise(
    n_total = n(),
    n_outliers = sum(is_outlier),
    ratio = n_outliers / n_total
  )


#What are these outliers in terms of topical study? 

bad_ratio <- bad_samples_info %>%
  count(project_topical_study_label, name = "n_bad") %>%
  left_join(
    pca_input_surface %>%
      count(project_topical_study_label, name = "n_total"),
    by = "project_topical_study_label"
  ) %>%
  mutate(ratio = n_bad / n_total) %>%
  arrange(desc(ratio))

bad_ratio


#Figure for th epaper with magma and distance to the centroid:
pcoa_df2 <- pcoa_df2 %>%
  left_join(intra_dist, by = "sample_id")

library(viridis)
p <- ggplot(pcoa_df2, aes(
  PCoA1, PCoA2,
  color = intra
)) +
  geom_segment(
    aes(xend = cx, yend = cy, group = sample_id),
    color = "grey70",
    alpha = 0.4,
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  geom_point(size = 2.5) +
  
  scale_color_viridis(
    option = "rocket",
    direction = -1,
    begin = 0,
    end = 0.85 ,
    name = "Distance to\ncentroid"
  ) +
  labs(
    x = paste0("PCoA1 (", pc1_var, "%)"),
    y = paste0("PCoA2 (", pc2_var, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 

p
pdf(
  file.path(
    output_path,
    "PCoA_Distance_to_centroid.pdf"
  ),
  width = 12,
  height = 5
)

print(p)

dev.off()

####3) Filtering #####

good_samples <- setdiff(
  unique(pca_input_surface$ESSENTIAL_sample_id),
  bad_samples$sample_id
)

pca_input_filtered <- pca_input_surface %>%
  filter(ESSENTIAL_sample_id %in% good_samples)

pca_matrix_filtered <- pca_input_filtered %>%
  select(starts_with("FT")) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

bray_dist_filtered <- vegdist(pca_matrix_filtered, method = "bray")

pcoa_res_filtered <- cmdscale(bray_dist_filtered, eig = TRUE, k = 2, add = TRUE)

pcoa_df_filtered <- as.data.frame(pcoa_res_filtered$points)
colnames(pcoa_df_filtered) <- c("PCoA1", "PCoA2")

pcoa_df_filtered$sample_id <- pca_input_filtered$ESSENTIAL_sample_id
pcoa_df_filtered$file_name <- pca_input_filtered$file_name
pcoa_df_filtered$station <- pca_input_filtered$ESSENTIAL_station
pcoa_df_filtered$topical_study <- pca_input_filtered$project_topical_study_label
pcoa_df_filtered$sample_id <- pca_input_filtered$ESSENTIAL_sample_id
pcoa_df_filtered$file_name <- pca_input_filtered$file_name
pcoa_df_filtered$station <- pca_input_filtered$ESSENTIAL_station
pcoa_df_filtered$topical_study <- pca_input_filtered$project_topical_study_label

centroids_filtered <- pcoa_df_filtered %>%
  group_by(sample_id) %>%
  summarise(
    cx = mean(PCoA1),
    cy = mean(PCoA2),
    .groups = "drop"
  )

pcoa_df2_filtered <- left_join(pcoa_df_filtered, centroids_filtered, by = "sample_id")

intra_dist_filtered <- pcoa_df2_filtered %>%
  group_by(sample_id) %>%
  summarise(
    intra = mean(sqrt((PCoA1 - cx)^2 + (PCoA2 - cy)^2))
  )

pcoa_df2_filtered <- pcoa_df2_filtered %>%
  left_join(intra_dist_filtered, by = "sample_id")


eig_vals_filt <- pcoa_res_filtered$eig
eig_vals_filt <- eig_vals_filt[eig_vals_filt > 0]

var_explained_filt <- eig_vals_filt / sum(eig_vals_filt)

pc1_var_filt <- round(var_explained_filt[1] * 100, 1)
pc2_var_filt <- round(var_explained_filt[2] * 100, 1)


p_filtered <- ggplot(pcoa_df2_filtered, aes(
  PCoA1, PCoA2
)) +
  geom_segment(
    aes(xend = cx, yend = cy, group = sample_id),
    color = "grey1",
    alpha = 0.4,
    linewidth = 0.3
  ) +
  geom_point(size = 2.5, 
             color = "grey35",
             alpha = 0.40) + 
  
  scale_color_viridis(
    option = "rocket",
    direction = -1,
    begin = 0,
    end = 0.85,
    name = "Distance to\ncentroid"
  ) +
  
  labs(
    x = paste0("PCoA1 (", pc1_var_filt, "%)"),
    y = paste0("PCoA2 (", pc2_var_filt, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 


p_filtered


####4) No filtering but average #####

pca_matrix_mean <- pca_input_surface %>%
  group_by(ESSENTIAL_sample_id) %>%
  summarise(across(starts_with("FT"), ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup()


metadata_mean <- pca_input_surface %>%
  group_by(ESSENTIAL_sample_id) %>%
  summarise(
    file_name = first(file_name),
    station = first(ESSENTIAL_station_label),
    topical_study = first(project_topical_study_label),
    .groups = "drop"
  )

pca_matrix_mean_only <- pca_matrix_mean %>%
  select(starts_with("FT")) %>%
  as.data.frame()

bray_dist_mean <- vegdist(pca_matrix_mean_only, method = "bray")

pcoa_res_mean <- cmdscale(bray_dist_mean, eig = TRUE, k = 2, add = TRUE)

pcoa_df_mean <- as.data.frame(pcoa_res_mean$points)
colnames(pcoa_df_mean) <- c("PCoA1", "PCoA2")


pcoa_df_mean$sample_id <- pca_matrix_mean$ESSENTIAL_sample_id

pcoa_df_mean <- left_join(
  pcoa_df_mean,
  metadata_mean,
  by = c("sample_id" = "ESSENTIAL_sample_id")
)



# On récupère les coordonnées de référence (non filtré)
ref <- pcoa_df %>%
  group_by(sample_id) %>%
  summarise(
    ref_PC1 = mean(PCoA1),
    ref_PC2 = mean(PCoA2),
    .groups = "drop"
  )

# Merge avec la version mean
align_df <- pcoa_df_mean %>%
  left_join(ref, by = "sample_id")

# Vérifier le signe via corrélation
cor_pc1 <- cor(align_df$PCoA1, align_df$ref_PC1, use = "complete.obs")
cor_pc2 <- cor(align_df$PCoA2, align_df$ref_PC2, use = "complete.obs")

# Flip si nécessaire
if (cor_pc1 < 0) {
  pcoa_df_mean$PCoA1 <- -pcoa_df_mean$PCoA1
}

if (cor_pc2 < 0) {
  pcoa_df_mean$PCoA2 <- -pcoa_df_mean$PCoA2
}


#pcoa_df_mean


eig_vals_mean <- pcoa_res_mean$eig
eig_vals_mean <- eig_vals_mean[eig_vals_mean > 0]

var_explained_mean <- eig_vals_mean / sum(eig_vals_mean)

pc1_var_mean <- round(var_explained_mean[1] * 100, 1)
pc2_var_mean <- round(var_explained_mean[2] * 100, 1)


p_mean <- ggplot(pcoa_df_mean, aes(
  PCoA1, PCoA2
)) +
  geom_point(size = 2.5, 
             color = "grey35",
             alpha = 0.40) + 
  
  scale_color_viridis(
    option = "rocket",
    direction = -1,
    begin = 0,
    end = 0.85,
    name = "Distance to\ncentroid"
  ) +
  
  labs(
    x = paste0("PCoA1 (", pc1_var_mean, "%)"),
    y = paste0("PCoA2 (", pc2_var_mean, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 

p_mean


library(patchwork)

p | p_filtered | p_mean

combined_plot <- (p | p_filtered | p_mean) +
  plot_annotation(tag_levels = "a")
ggsave(
  filename = file.path(
    output_path,
    "Figure_SX_PCoA_processing.pdf"
  ),
  plot = combined_plot,
  width = 12,
  height = 5,
  device = cairo_pdf
)

####5)Effect MeOH####



colnames(metadata)

#I want to compare Evaporation_status with Code_inlab_operator and ESSENTIAL_event_latitude_N_midpoint

metadata_mean <- pca_input_surface %>%
  group_by(ESSENTIAL_sample_id) %>%
  summarise(
    file_name = first(file_name),
    station = first(ESSENTIAL_station_label),
    topical_study = first(project_topical_study_label),
    Evaporation_status = first(Evaporation_status),
    Code_inlab_operator = first(Code_inlab_operator),
    ESSENTIAL_event_latitude_N_midpoint = first(ESSENTIAL_event_latitude_N_midpoint),
    .groups = "drop"
  )

metadata_mean %>%
  count(Code_inlab_operator, Evaporation_status)

metadata$event_label
pcoa_df_mean <- left_join(
  pcoa_df_mean,
  metadata_mean,
  by = c("sample_id" = "ESSENTIAL_sample_id")
)

unique(pcoa_df_mean$Evaporation_status)
p_mean_evaporation <- ggplot(pcoa_df_mean, aes(
  PCoA1, PCoA2,
  color = Evaporation_status  
)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_viridis_d(
    option = "mako",
    direction = -1,
    begin = 0,
    end = 0.90,
    name = "Evaporation\nstatus"
  ) +
  
  labs(
    x = paste0("PCoA1 (", pc1_var_mean, "%)"),
    y = paste0("PCoA2 (", pc2_var_mean, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 

p_mean_evaporation

p_mean_operator <- ggplot(pcoa_df_mean, aes(
  PCoA1, PCoA2,
  color = Code_inlab_operator  
)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_viridis_d(
    option = "mako",
    direction = -1,
    begin = 0.4,
    end = 0.90,
    name = "Operator"
  ) +
  
  
  labs(
    x = paste0("PCoA1 (", pc1_var_mean, "%)"),
    y = paste0("PCoA2 (", pc2_var_mean, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 


p_mean_operator


p_mean_latitude <- ggplot(pcoa_df_mean, aes(
  PCoA1, PCoA2,
  color = ESSENTIAL_event_latitude_N_midpoint 
)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_viridis_c(
    option = "mako",
    direction = 1,
    begin = 0,
    end = 0.90,
    name = "Latitude"
  ) +
  
  labs(
    x = paste0("PCoA1 (", pc1_var_mean, "%)"),
    y = paste0("PCoA2 (", pc2_var_mean, "%)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold")
  ) 

p_mean_latitude


p_mean_evaporation | p_mean_operator | p_mean_latitude

combined_plot <- (p_mean_evaporation | p_mean_operator | p_mean_latitude) +
  plot_annotation(tag_levels = "a")

ggsave(
  filename = file.path(
    output_path,"Figure_SX_PCoA_covariants.pdf"),
  plot = combined_plot,
  width = 12,
  height = 5,
  device = cairo_pdf
)


cor(metadata$ESSENTIAL_event_latitude_N_midpoint,
    as.numeric(metadata$Evaporation_status),
    method = "spearman")

#monotonic relashionship?7
ggplot(metadata, aes(
  x = ESSENTIAL_event_latitude_N_midpoint,
  y = as.numeric(Evaporation_status)
)) +
  geom_jitter(height = 0.1, width = 0) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_classic()





