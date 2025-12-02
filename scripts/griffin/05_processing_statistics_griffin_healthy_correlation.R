library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)
library(readxl)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(gridExtra)
library(data.table)
library(psych)
library(here)


output_dir <- here::here("data", "griffin", "01_processing_statistics_griffin_scores")
output_dir_combined <- here::here("data", "griffin")
date <- Sys.Date();

# TSV files option
samples_path <- here::here("raw_data", "cohort", "HBOC_cohort_metadata.xlsx") 
data_samples <- read_excel(samples_path, sheet="April22")

data_samples <- data_samples[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_samples$`Library Name`), ]
data_samples <- data_samples[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_samples$`Library Name`), ]

# R Objects option - since data is large
path <- here::here("raw_data", "griffin", "nucleosome_accessibility_HBOC.rds")

# robject - since data is large
healthy_path <- here::here("raw_data", "griffin", "nucleosome_accessibility_HBC.rds")


# read in robjects 
data_V1 <- readRDS(path)
data_V1 <- data_V1[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_V1$Sample), ]
data_V1 <- data_V1[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_V1$Sample), ]

healthy_data_V1 <- readRDS(healthy_path)

# Define the function to calculate medians for rest of cohort. 
calculate_median_sd <- function(category, label, data_samples, data_griffin) {
  # Filter samples
  samples <- data_samples[data_samples$source == category, ]

  # Select relevant columns in data_griffin
  griffin_data <- data_griffin[, c("distance", samples$'Library Name')]
  
  # Calculate medians and standard deviations
  medians <- rowMedians(as.matrix(griffin_data[, !(colnames(griffin_data) == "distance")]))
  sds <- rowSds(as.matrix(griffin_data[, !(colnames(griffin_data) == "distance")]))
  
  # Create the summary dataframe
  summary_df <- as.data.frame(cbind(distance, medians, sds))
  colnames(summary_df) <- c("distance", "median", "sd")
  summary_df$diag <- label

  
  return(summary_df)
}

sites <- c("LGG", "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", 
           "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO", "PCPG", "PRAD", "SKCM", "STAD", 
           "TGCT", "THCA", "UCEC", "OV")

names <- c("Brain Lower Grade Glioma", "Adrenocortical Carcinoma", "Bladder Cancer", "Breast Cancer", 
           "Cervical and Endocervical Cancer", "Cholangiocarcinoma", 
           "Colon Cancer", "Esophageal Cancer", "Glioblastoma", "Head and Neck Cancer", 
           "Kidney Clear Cell Carcinoma", "Kidney Papillary Cell Carcinoma", "Liver Cancer",
           "Lung Adenocarcinoma", "Lung Squamous Cell Carcinoma", "Mesothelioma", 
           "Pheochromocytoma and Paraganglioma", "Prostate Cancer", "Melanoma", 
           "Stomach Cancer", "Testicular Germ Cell Tumors", "Thyroid Cancer", "Uterine Endometrial Cancer", "Ovarian Cancer")

lists <- c("brain", "adrenocortical", "bladder", "breast", "cervical_endocervical", "cholangiocarcinoma", 
           "colon", "esophageal", "glioblastoma", "head_neck",
           "kidney_clear_cell", "kidney_papillary_cell", "liver", "lung_adenocarcinoma", 
           "lung_squamous_cell", "mesothelioma", "pheochromocytoma_paraganglioma", "prostate", 
           "melanoma", "stomach", "testicular_germ_cell", "thyroid", "uterine_endometrial", "ovarian")

# Initialize empty lists to store dataframes for cancer and normal samples
all_cancer_list <- list()
all_normal_list <- list()

for (i in c(1:length(sites))) {
  
  data <- data_V1
  healthy_data <- healthy_data_V1
  ### Set variables
  site <- sites[[i]]
  name <- names[[i]]
  list <- lists[[i]]
  
  head(data)
  # Ensure 'site_name' column exists before filtering
  if (!("site_name" %in% colnames(data)) || !("site_name" %in% colnames(healthy_data))) {
    stop("Column 'site_name' not found in data or healthy_data")
  }
  
  # Alternative to filter function using base R subsetting
  site_data <- data[data$site_name == site, ]
  site_healthy_data <- healthy_data[healthy_data$site_name == site, ]
  
  
  # Convert from long to wide format
  site_data_wide <- site_data %>%
    pivot_wider(names_from = Sample, values_from = Coverage)
  
  site_healthy_data_wide <- site_healthy_data %>%
    pivot_wider(names_from = Sample, values_from = Coverage)
  
  # Remove unnecessary columns by explicitly dropping them
  data_griffin <- site_data_wide[, !(colnames(site_data_wide) %in% c("site_name", "site_type"))]
  data_normal <- site_healthy_data_wide[, !(colnames(site_healthy_data_wide) %in% c("site_name", "site_type"))]
  
  # Rename Position to distance
  colnames(data_griffin)[colnames(data_griffin) == "Position"] <- "distance"
  colnames(data_normal)[colnames(data_normal) == "Position"] <- "distance"
  
  ### Apply Savitzky-Golay filter
  distance <- data_griffin$distance
  data_griffin <- sapply(data_griffin[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_griffin <- as.data.frame(data_griffin)
  
  data_normal <- sapply(data_normal[, -1], function(x) savgol(x, 11, forder = 3, dorder = 0))
  data_normal <- as.data.frame(data_normal)
  
  ### Adjust and center curves
  means <- colMeans(data_griffin[c(54,55,56,77,78,79), ])
  means <- means - mean(means)
  data_griffin <- sweep(data_griffin, 2, means)
  data_griffin$distance <- distance
  
  means <- colMeans(data_normal[c(54,55,56,77,78,79), ])
  means <- means - mean(means)
  data_normal <- sweep(data_normal, 2, means)
  data_normal$distance <- distance
  
  ### Order based on clinical information
  data_samples <- data_samples[data_samples$'Library Name' %in% colnames(data_griffin), ]
  data_samples <- data_samples[order(data_samples$timepoint), ]
  data_samples$diag <- "HBOC"
  ### Format dataframes
  data_griffin <- data_griffin[, c("distance", data_samples$'Library Name')]
  
  # Filter data between -50 and 50 for both data_normal and data_griffin
  data_normal_filtered <- data_normal[data_normal$distance >= -30 & data_normal$distance <= 30, ]
  data_griffin_filtered <- data_griffin[data_griffin$distance >= -30 & data_griffin$distance <= 30, ]
  
  # Extract 'distance' column for later use
  distance <- data_normal_filtered$distance
  
  ### Calculate the healthy median
  normal_median <- rowMedians(as.matrix(data_normal_filtered[, !(colnames(data_normal_filtered) == "distance")]))
  normal_sd <- rowSds(as.matrix(data_normal_filtered[, !(colnames(data_normal_filtered) == "distance")]))
  
  ### Make median tables - Healthy
  normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
  colnames(normal_median) <- c("distance", "median", "sd")
  normal_median$diag <- "Healthy"
  
  # 1. get mean differences 
  mean_diff_scores <- colMeans(sweep(data_griffin_filtered[, !(colnames(data_griffin_filtered) == "distance")], 1, normal_median$median, FUN = "-"))
  # 2. Calculate mean differences for normal patients using median
  normal_mean_diff <- colMeans(sweep(data_normal_filtered[,!(colnames(data_normal_filtered) == "distance")], 1, normal_median$median, FUN = "-"))
  # 3. Determine the 95th percentile as the cutoff
  cutoff <- quantile(normal_mean_diff, 0.95)
  # 4. Check which cancer patients exceed this cutoff
  significant_patients_mean_diff <- mean_diff_scores > cutoff
  
  
  # Remove the distance column from both cancer and normal data
  data_griffin_filtered_no_dist <- data_griffin_filtered[, !names(data_griffin_filtered) %in% "distance"]
  data_normal_filtered_no_dist <- data_normal_filtered[, !names(data_normal_filtered) %in% "distance"]
  # 1. Calculate the mean and standard deviation for each site (row) using the normal samples
  normal_means <- apply(data_normal_filtered_no_dist, 1, mean)
  normal_sds <- apply(data_normal_filtered_no_dist, 1, sd)
  # 2. Calculate z-scores for each cancer patient
  z_scores <- sweep(data_griffin_filtered_no_dist, 1, normal_means, FUN = "-")
  z_scores <- sweep(z_scores, 1, normal_sds, FUN = "/")
  
  
  z_scores_normal <- sweep(data_normal_filtered_no_dist, 1, normal_means, FUN = "-")
  z_scores_normal <- sweep(z_scores_normal, 1, normal_sds, FUN = "/")
  
  # 3. Calculate the average z-score for each cancer patient across the 5 sites
  average_z_scores <- colMeans(z_scores)
  average_z_scores_normal <- colMeans(z_scores_normal, na.rm = TRUE)
  
  
  # Set the significance threshold for z-scores
  threshold <- -1.96
  # Identify patients with significant z-scores (absolute value greater than 1.96)
  significant_patients_zscores <- average_z_scores[average_z_scores < threshold]
  cat(list)
  print(length(significant_patients_zscores))

  # Create the cancer patients' data frame with mean differences and z-scores
  cancer_df <- data.frame(
    Patient = names(mean_diff_scores),
    MeanDifference = mean_diff_scores,
    AverageZScore = average_z_scores
  )
  colnames(cancer_df)[2:3] <- paste0(c("MeanDifference_", "AverageZScore_"), list)
  # Create the normal patients' data frame with normal mean differences
  normal_df <- data.frame(
    Patient = names(normal_mean_diff),
    NormalMeanDifference = normal_mean_diff, 
    NormalAverageZScore = average_z_scores_normal
  )
  colnames(normal_df)[2:3] <- paste0(c("MeanDifference_", "AverageZScore_"), list)
  
  # Collect results for all sites
  all_cancer_list[[i]] <- cancer_df
  all_normal_list[[i]] <- normal_df
  
  # Print completion message for debugging
  cat("Completed processing for site:", list, "\n")
  
} 

# Combine cancer dataframes by "Patient" column (wide format)
all_cancer_df <- Reduce(function(x, y) merge(x, y, by = "Patient", all = TRUE), all_cancer_list)

# Combine normal dataframes by "Patient" column (wide format)
all_normal_df <- Reduce(function(x, y) merge(x, y, by = "Patient", all = TRUE), all_normal_list)

# Write combined CSV files
write.csv(all_cancer_df, file = file.path(output_dir_combined, paste0("05_processing_griffin_nucleosome_accesibility_all_cancer_scores", ".csv")), row.names = FALSE)
write.csv(all_normal_df, file = file.path(output_dir_combined, paste0("05_processing_griffin_nucleosome_accesibility_all_normal_scores", ".csv")), row.names = FALSE)
  