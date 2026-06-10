#I used MetabCombiner to link together the full MS and the DDA dataset
#I followed the instructions provided in this vignette:
#https://www.bioconductor.org/packages/release/bioc/vignettes/metabCombiner/inst/doc/metabCombiner_vignette.html
#Please consult it for more info!

####1) Installation ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


#Prepare your environment
library(here)
base_path <- here()
input_path <- file.path(base_path, "1_Input_files")
plot_path <- file.path(base_path, "2_Quality_graphs")
out_path <- file.path(base_path, "3_Output_files")
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

#BiocManager::install("metabCombiner") #if necessary

#Loading of the packages:
library(metabCombiner)
library(RColorBrewer)
library(pander)
library(magrittr)
library(SummarizedExperiment)
library(ggplot2)
library(tidyverse)
library(readxl)
library(tidyr)
library(dplyr)

# feature table full MS
feature_table_fullMS <- read.table(
  file.path(input_path,
            "feature_table_MS1_positive1.txt"),
  header = TRUE
)
colnames(feature_table_fullMS)
colnames(feature_table_fullMS)[which(names(feature_table_fullMS) == "Row.names")] <- "featureId"
head(colnames(feature_table_fullMS))
feature_table_fullMS <- feature_table_fullMS %>%
  select(featureId, mzmed, rtmed, starts_with("X"))



# feature table DDA
feature_table_DDA <- read.table(
  file.path(input_path,
            "DDAfeature_table_positive1.txt"),
  header = TRUE
)
colnames(feature_table_DDA)
colnames(feature_table_DDA)[which(names(feature_table_DDA) == "Row.names")] <- "featureId"
feature_table_DDA <- feature_table_DDA %>%
  select(featureId, mzmed, rtmed, starts_with("X"))


####2) Data formatting and filtering ####


#####2.1) Feature table full MS ####

feature_table_fullMS$rtmed <- feature_table_fullMS$rtmed / 60 #RT should be in minutes
head(sort(feature_table_fullMS$rtmed), 10) 
tail(sort(feature_table_fullMS$rtmed), 10) 
head(colnames(feature_table_fullMS))
pfullMS <- metabData(table = feature_table_fullMS, 
                 mz = "mzmed", 
                 rt = "rtmed", 
                 id = "featureId", 
                 adduct = NULL, 
                 samples = "^X", #Column names starting with X
                 extra = NULL,  
                 rtmin = "min", 
                 rtmax = "max", #In min!
                 misspc = 99.9,
                 measure = "median", 
                 zero = FALSE, 
                 duplicate = c(0.0025, 0.05))
getSamples(pfullMS) ##should print out my samples names
getStats(pfullMS) ##prints a list of dataset statistics:
#input_size:99207
#filtered_by_rt:0
#filtered_as_duplicates: 197
#filtered_by_missingness: 196
#final_count: 98814
print(pfullMS)   ##object summary
#Total Samples: 1393   Total Extra: 0 
#Mass Range: 100.011-1492.044 Da
#Time Range: 1.0017-11.6257 minutes
#Input Feature Count: 99207 
#Final Feature Count: 98814 

#####2.2) Feature table DDA ####

feature_table_DDA$rtmed <- feature_table_DDA$rtmed / 60
head(sort(feature_table_DDA$rtmed), 10) #Normal
tail(sort(feature_table_DDA$rtmed), 10) #Fall to far apart, should stop at 11.54
pDDA <- metabData(table = feature_table_DDA, 
                     mz = "mzmed", 
                     rt = "rtmed", 
                     id = "featureId", 
                     adduct = NULL, 
                     samples = "^X", #Column names starting with X
                     extra = NULL,  
                     rtmin = "min", 
                     rtmax = 11.54, #In min!
                     misspc = 99.9,
                     measure = "median", 
                     zero = FALSE, 
                     duplicate = c(0.0025, 0.05))
getStats(pDDA) ##prints a list of dataset statistics:
#input_size:114571
#filtered_by_rt:3
#filtered_as_duplicates: 327
#filtered_by_missingness: 0
#final_count: 114241
print(pDDA)   ##object summary
#Total Samples: 411   Total Extra: 0 
#Mass Range: 80.0343-1498.798 Da
#Time Range: 0.9879-11.5317 minutes
#Input Feature Count: 114571
#Final Feature Count: 114241

####3) Feature m/z Grouping and Pairwise Alignment Detection ####

#pFullMS possesses the longest RT range -> X
#pDDA possesses the shortest RT range -> Y

p.combined = metabCombiner(xdata = pfullMS, ydata = pDDA, binGap = 0.005)
p.results = combinedTable(p.combined)

####4) Anchor Selection and RT Mapping ####

#####4.1) Anchor Selection ####

p.combined.2 = selectAnchors(p.combined, 
                             windx = 0.03,
                             windy = 0.02, 
                             tolQ = 0.3, 
                             tolmz = 0.003, 
                             tolrtq = 0.3, 
                             useID = FALSE)

a = getAnchors(p.combined.2)
plot(a$rtx, a$rty, main = "Fit Template", xlab = "rtx", ylab = "rty")

p <- ggplot(a, aes(x = rtx, y = rty)) +
  geom_point(
    color = "black",
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    title = "Retention Time Alignment Anchors",
    x = "Retention time in min (Full MS)",
    y = "Retention time in min (DDA)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 12, by = 1),       # major ticks
    minor_breaks = seq(0, 12, by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 12, by = 1),
    minor_breaks = seq(0, 12, by = 0.5)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )


ggsave(
  filename = file.path(
    plot_path,
    "RT_alignment_anchors_metabCombiner.pdf"
  ),
  plot = p,
  width = 7,
  height = 6,
  device = cairo_pdf
)
#####4.2) Model-fitting####

set.seed(100) #controls cross validation pseudo-randomness
p.combined.3 = fit_gam(p.combined.2, 
                       useID = FALSE, 
                       k = seq(12,20,2), 
                       iterFilter = 2, 
                       coef = 2, 
                       prop = 0.5, 
                       bs = "bs", 
                       family = "scat", 
                       m = c(3,2))
pdf(
  file.path(plot_path,
            "RT_projection_model_clean.pdf"),
  width = 7,
  height = 6
)

par(
  bty = "l",
  cex.lab = 1.3,
  cex.axis = 1.2,
  font.lab = 2
)

plot(
  p.combined.3,
  main = "GAM Retention Time Mapping",
  xlab = "Retention time (Full MS) [min]",
  ylab = "Retention time (DDA) [min]",
  lcol = "#2166AC",
  pcol = "black",
  lwd = 3,
  pch = 16,
  outlier = "highlight"
)

dev.off()


#####4.3) Feature Pair Alignment Scoring####

p.combined.4 <- calcScores(p.combined.3,
                           A = 0.03,
                           B = 15,
                           C = 0.5,
                           usePPM = TRUE,
                           useAdduct = FALSE)


#Clean your environment to delete everything except p.combined.4

####5) Table Annotation & Reduction####
save.image(
  file = file.path(
    out_path,
    "EnvironmentLight.RData"
  )
)

combined.table = combinedTable(p.combined.4)
rm(p.combined.4)
gc()
write.table(combined.table,
            file = file.path(out_path,
                             "Combined.Table.Report.txt"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
gc()

#Score-based conflict detection
combined.table.2 = labelRows(combined.table, 
                             minScore = 0.5, 
                             maxRankX = 3, 
                             maxRankY = 3, 
                             method = "score", 
                             delta = 0.2, 
                             remove = FALSE, 
                             balanced = TRUE)

getStats(combined.table.2)
write.table(
  combined.table.2,
  file = file.path(
    out_path,
    "CombinedTable_withQuality.txt"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)


combined.table.3 <- combined.table.2[, !grepl("^X", colnames(combined.table.2))]
table(combined.table.3$labels)
write.table(combined.table.3,
            file = file.path(out_path,
                             "CombinedTable_withQuality_Short.txt"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")


combined.table.4 <- combined.table.3[!combined.table.3$labels %in% c("REMOVE", "CONFLICT"), ]
write.table(combined.table.4,
            file = file.path(out_path,
                             "CombinedTable_withQuality_Minimal.txt"),
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
