library(gsignal)
library(dplyr)
library(matrixStats)
library(cowplot)
library(GeneCycle)
library(pracma)
library(readxl)
library(tidyr)
library(ggplot2)
library(here)

# Project root is automatically detected by here()
cat("Project root:", here::here(), "\n")


# Output directory
outdir <- here::here("figures", "griffin", "02_figure_griffin_by_source")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Metadata path
samples_path <- here::here("raw_data", "cohort", "HBOC_cohort_metadata.xlsx")

data_samples <- read_excel(samples_path, sheet = "April22")
data_samples <- data_samples[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_samples$`Library Name`), ]
data_samples <- data_samples[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_samples$`Library Name`), ]

# R Objects option - since data is large
path <- here::here("raw_data", "griffin", "nucleosome_accessibility_HBOC.rds")

healthy_path <- here::here("raw_data", "griffin", "nucleosome_accessibility_HBC.rds")

# Read in R objects
data_V1 <- readRDS(path)
data_V1 <- data_V1[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_V1$Sample), ]
data_V1 <- data_V1[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_V1$Sample), ]

healthy_data_V1 <- readRDS(healthy_path)

# Define the function to calculate medians for rest of cohort. 
calculate_median_sd <- function(category, label, data_samples, data_griffin) {
  # Filter samples
  samples <- data_samples[data_samples$source == category, ]
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
sites <- c("OV", "BRCA", "PRAD", "LUAD", "BLCA", "LGG", "ACC", "CESC", "CHOL", "COAD", "ESCA", "GBM", "HNSC", 
           "KIRC", "KIRP", "LIHC", "LUSC", "MESO", "PCPG", "SKCM", "STAD", "TGCT", "THCA", "UCEC")

names <- c("Ovarian Cancer", "Breast Cancer", "Prostate Cancer", "Lung Adenocarcinoma", "Bladder Cancer", 
           "Brain Lower Grade Glioma", "Adrenocortical Carcinoma", "Cervical and Endocervical Cancer", 
           "Cholangiocarcinoma", "Colon Cancer", "Esophageal Cancer", "Glioblastoma", "Head and Neck Cancer", 
           "Kidney Clear Cell Carcinoma", "Kidney Papillary Cell Carcinoma", "Liver Cancer", 
           "Lung Squamous Cell Carcinoma", "Mesothelioma", "Pheochromocytoma and Paraganglioma", "Melanoma", 
           "Stomach Cancer", "Testicular Germ Cell Tumors", "Thyroid Cancer", "Uterine Endometrial Cancer")

lists <- c("ovarian", "breast", "prostate", "lung_adenocarcinoma", "bladder", "brain", "adrenocortical", 
           "cervical_endocervical", "cholangiocarcinoma", "colon", "esophageal", "glioblastoma", "head_neck", 
           "kidney_clear_cell", "kidney_papillary_cell", "liver", "lung_squamous_cell", "mesothelioma", 
           "pheochromocytoma_paraganglioma", "melanoma", "stomach", "testicular_germ_cell", "thyroid", "uterine_endometrial")


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
  
  ### Calculate the healthy median
  normal_median <- rowMedians(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  normal_sd <- rowSds(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
  
  ### Make median tables - Healthy 
  normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
  colnames(normal_median) <- c("distance", "median", "sd")
  normal_median$diag <- "Healthy Control"
  
  # Call the function for each category
  cancer_status_negative_survivor_no_median <- calculate_median_sd("cancer_status_negative_survivor_no", "Cancer Negative\n (Survivor No)", data_samples, data_griffin)
  cancer_status_negative_survivor_yes_median <- calculate_median_sd("cancer_status_negative_survivor_yes", "Cancer Negative\n (Survivor)", data_samples, data_griffin)
  first_positive_survivor_no_median <- calculate_median_sd("first_positive_survivor_no", "Cancer Positive\n (First Positive)", data_samples, data_griffin)
  cancer_status_positive_survivor_yes_median <- calculate_median_sd("cancer_status_positive_survivor_yes", "Cancer Positive\n (Survivor)", data_samples, data_griffin)
  
  # bind all the dataframes. 
  data_median <- bind_rows(cancer_status_positive_survivor_yes_median, first_positive_survivor_no_median, cancer_status_negative_survivor_yes_median, cancer_status_negative_survivor_no_median, normal_median)
  
  ## ++++++++++ Important step +++++++++++++++++++
  samples <- data_samples
  ## ++++++++++ Important step +++++++++++++++++++
  
  data <- data_griffin[, c("distance", samples$'Library Name')]
  #data$distance <- data$distance - 8.5
  data_melt <- reshape2::melt(data, id = "distance")
  data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "Library Name", all = TRUE)
  data_melt <- data_melt[data_melt$distance > -901 & data_melt$distance < 901, ]
  
  # Use factor to set levels and labels for a new column diag
  data_melt$diag <- factor(data_melt$source, 
                           levels = c("first_positive_survivor_no", "cancer_status_positive_survivor_yes", "cancer_status_negative_survivor_no", "cancer_status_negative_survivor_yes"),
                           labels = c("Cancer Positive\n (First Positive)", "Cancer Positive\n (Survivor)", "Cancer Negative\n (Survivor No)", "Cancer Negative\n (Survivor)"))
  
  ### Set median line
  median <- min(data_median$median[data_median$diag == "Healthy Control" & data_median$distance %in% c("-30", "-15", "0", "15", "30")])
  
  ### Set Theme
  theme <- theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
                      axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      legend.position = "right",
                      legend.key = element_rect(fill = "white", colour = NA),  # remove legend key border
                      legend.text = element_text(size = 10),
                      legend.title = element_text(size = 12),
                      strip.background = element_rect(fill = "white"),
                      strip.text = element_text(size = 13),
                      axis.text = element_text(size = 8),
                      axis.title = element_text(size = 13),
                      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),  
                      axis.text.y = element_text(size = 12)
                 )
  
  ### Plot curves
  fig <- ggplot(data_median) +
    geom_line(aes(distance, median), color = "black", size = 0.5) +
    geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd), fill = "black", alpha = 0.5) +
    geom_line(data = data_melt, aes(distance, value, color=diag, group = variable), alpha = 0.5, size = 0.5) +
    geom_hline(yintercept = median, linetype = "dashed", size = 0.5) +
    ggtitle(paste("Nucleosome Accesbility Coverage by Cancer Group for", name)) +
    xlab("Distance from site (bp)") + 
    ylab("Coverage") +
    labs(color = "Cancer Groups", fill = "") +
    guides(alpha = "none", color = guide_legend(override.aes = list(linewidth =3))) +
    scale_color_manual(labels = c("cancer_status_negative_survivor_yes"= "Cancer Negative\n (Survivor)", 
                                  "cancer_status_positive_survivor_yes" = "Cancer Positive\n (Survivor)", 
                                  "first_positive_survivor_no" = "Cancer Positive\n (First Positive)",
                                  "cancer_status_negative_survivor_no"= "Cancer Negative\n (Survivor No)"),
                       values = c("Cancer Negative\n (Survivor)"="#1f77b4","Cancer Positive\n (Survivor)" = "#f4c842",
                                  "Cancer Positive\n (First Positive)" = "#f87544",
                                  "Cancer Negative\n (Survivor No)"="purple")) +
    #scale_alpha_manual(values = c("negative" = 0.25, "postiive" = 0.75)) +
    facet_wrap(.~diag, nrow=1) +
    theme +
    scale_x_continuous(limits = c(-900, 900), expand = c(0,0))
  fig
  
  ggsave(file.path(outdir, paste0("02_griffin_nucleosome_curves_TCGA_by_source_", site, ".pdf")), fig, width = 14, height = 6)
  plot_name <- paste0(site, "_plot")
  assign(plot_name, fig)
  
  # Print completion message for debugging
  cat("Completed processing for site:", site, "\n")
  
} 
