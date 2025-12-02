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

# ==========================
# Output Directory
# ==========================
outdir <- here::here("figures", "griffin", "03_figure_griffin_by_cancer")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================
# Load Sample Metadata
# ==========================
samples_path <- here::here("raw_data", "cohort", "HBOC_cohort_metadata.xlsx")
data_samples <- read_excel(samples_path, sheet = "April22")

# Remove failed samples
data_samples <- data_samples[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_samples$`Library Name`), ]
data_samples <- data_samples[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_samples$`Library Name`), ]

# Standardize cancer_type
data_samples$cancer_type[data_samples$cancer_type == "hgsoc"] <- "ovarian"
data_samples$cancer_type[data_samples$cancer_type == "breast, hgsoc"] <- "ovarian"
data_samples$cancer_type[data_samples$cancer_type == "thyroid, hgsoc"] <- "ovarian"

# ==========================
# Load RDS Data
# ==========================
path <- here::here("raw_data", "griffin","nucleosome_accessibility_HBOC.rds")

healthy_path <- here::here("raw_data", "griffin","nucleosome_accessibility_HBC.rds")

# Read data objects
data_V1 <- readRDS(path)
data_V1 <- data_V1[!grepl("CHARMQ_0788_Pl_T_PG_T-788", data_V1$Sample), ]
data_V1 <- data_V1[!grepl("CHARMQ_0207_Pl_T_PG_T-207", data_V1$Sample), ]

healthy_data_V1 <- readRDS(healthy_path)

# Prepare placeholders for the two target plots
UEC_breast_plot  <- NULL
UEC_ovarian_plot <- NULL

# Define the function to calculate medians for rest of cohort. 
calculate_median_sd <- function(category, label, data_samples, data_griffin) {
  # Filter samples
  samples <- data_samples[!is.na(data_samples$cancer_type) & data_samples$cancer_type == category, ]
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

cancer_type_patients <- c("ovarian","breast", "liver", "thyroid", "prostate")

for (t in seq_along(cancer_type_patients)) { 
  for (i in seq_along(sites)) {
    
    data <- data_V1
    healthy_data <- healthy_data_V1
    
    ### Set variables
    site <- sites[[i]]
    name <- names[[i]]
    list <- lists[[i]]
    cancer_type <- cancer_type_patients[[t]]
    pos <- paste0("BRCA1/2m Cancer Positive\n (Excluding ", cancer_type, ")")
    
    head(data)
    # Ensure 'site_name' column exists before filtering
    if (!("site_name" %in% colnames(data)) || !("site_name" %in% colnames(healthy_data))) {
      stop("Column 'site_name' not found in data or healthy_data")
    }
    
    # Subset by site
    site_data <- data[data$site_name == site, ]
    site_healthy_data <- healthy_data[healthy_data$site_name == site, ]
    

    
    # Convert from long to wide format
    site_data_wide <- site_data %>%
      pivot_wider(names_from = Sample, values_from = Coverage)
    
    site_healthy_data_wide <- site_healthy_data %>%
      pivot_wider(names_from = Sample, values_from = Coverage)
    
    # Remove unnecessary columns
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
    data_samples <- data_samples[data_samples$`Library Name` %in% colnames(data_griffin), ]
    data_samples <- data_samples[order(data_samples$timepoint), ]
    data_samples$diag <- ""
    ### Format dataframes
    data_griffin <- data_griffin[, c("distance", data_samples$`Library Name`)]
    
    ### Calculate the healthy median
    normal_median <- rowMedians(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
    normal_sd <- rowSds(as.matrix(data_normal[, !(colnames(data_normal) == "distance")]))
    
    ### Make median tables - Healthy 
    normal_median <- as.data.frame(cbind(distance, normal_median, normal_sd))
    colnames(normal_median) <- c("distance", "median", "sd")
    normal_median$diag <- "Healthy Control"
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ### Calculate the Negative median 
    # Filter for negative samples 
    samples_neg <- data_samples[!is.na(data_samples$cancer_status) & data_samples$cancer_status == "negative", ]
    # Select relevant columns in data_griffin
    griffin_data_neg <- data_griffin[, c("distance", data_samples$`Library Name`)]
    
    neg_median <- rowMedians(as.matrix(griffin_data_neg[, !(colnames(griffin_data_neg) == "distance")]))
    neg_sd <- rowSds(as.matrix(griffin_data_neg[, !(colnames(griffin_data_neg) == "distance")]))
    ### Make median tables - Negative 
    neg_median <- as.data.frame(cbind(distance, neg_median, neg_sd))
    colnames(neg_median) <- c("distance", "median", "sd")
    neg_median$diag <- "BRCA1/2m Cancer Negative"
    
    ### Calculate the positive (exluding specified cancer_type) median 
    samples_pos <- data_samples[!is.na(data_samples$cancer_type) & data_samples$cancer_status == "positive" & !data_samples$cancer_type == cancer_type, ]
    griffin_data_pos <- data_griffin[, c("distance", data_samples$`Library Name`)]
    
    pos_median <- rowMedians(as.matrix(griffin_data_pos[, !(colnames(griffin_data_pos) == "distance")]))
    pos_sd <- rowSds(as.matrix(griffin_data_pos[, !(colnames(griffin_data_pos) == "distance")]))
    ### Make median tables - Positive excluding specified cancer_type
    pos_median <- as.data.frame(cbind(distance, pos_median, pos_sd))
    colnames(pos_median) <- c("distance", "median", "sd")
    pos_median$diag <- pos
    
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Call the function for the specified cancer_type 
    cancer_type_median <- calculate_median_sd(cancer_type, cancer_type, data_samples, data_griffin)
    
    # bind all the dataframes. 
    data_median <- bind_rows(cancer_type_median, normal_median, neg_median, pos_median)
    
    ## ++++++++++ Important step +++++++++++++++++++
    samples <- data_samples
    ## ++++++++++ Important step +++++++++++++++++++
    
    data <- data_griffin[, c("distance", samples$`Library Name`)]
    
    # data$distance <- data$distance - 8.5
    data_melt <- reshape2::melt(data, id = "distance")
    data_melt <- merge(data_melt, samples, by.x = "variable", by.y = "Library Name", all = TRUE)
    data_melt <- data_melt[data_melt$distance > -901 & data_melt$distance < 901, ]
    
    # Set the diag column for different categories
    data_melt$diag[!is.na(data_melt$cancer_status) & data_melt$cancer_status == "negative"] <- "BRCA1/2m Cancer Negative"
    data_melt$diag[!is.na(data_melt$cancer_status) & data_melt$cancer_status == "positive"] <- pos
    data_melt$diag[!is.na(data_melt$cancer_type) & data_melt$cancer_type == cancer_type] <- cancer_type
    # Use factor to set levels and labels for a new column diag
    data_melt$diag <- factor(data_melt$diag,  
                             levels = c("BRCA1/2m Cancer Negative", cancer_type , pos, "Healthy Control"))
    
    # Drop rows with NA in the diag column
    data_melt <- data_melt[!is.na(data_melt$diag), ]
    
    ### Set median line
    median <- min(data_median$median[data_median$diag == "Healthy Control" & data_median$distance %in% c(-30, -15, 0, 15, 30)])
    
    ### Set Theme
    theme <- theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), 
                   axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   legend.position = "right",
                   legend.key = element_rect(fill = "white", colour = NA),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 12),
                   strip.background = element_rect(fill = "white"),
                   strip.text = element_text(size = 13),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 13),
                   axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
                   axis.text.y = element_text(size = 12))
    
    # Define color values dynamically
    color_values <- c("BRCA1/2m Cancer Negative" = "#76EE00",
                      setNames("#FFA500", cancer_type),
                      setNames("#EE0000", pos))
    ### Plot curves
    fig <- ggplot(data_median) +
      geom_line(aes(distance, median), color = "black", size = 0.5) +
      geom_ribbon(aes(distance, ymin = median - sd, ymax = median + sd), fill = "black", alpha = 0.5) +
      geom_line(data = data_melt, aes(distance, value, color=diag, group = variable), alpha = 0.5, size = 0.5) +
      geom_hline(yintercept = median, linetype = "dashed", size = 0.5) +
      ggtitle(paste("Nucleosome Accesbility Coverage by Cancer Type for", name)) +
      xlab("Distance from site (bp)") +
      ylab("Coverage") +
      labs(color = "Cancer Group", fill = "") +
      guides(alpha = "none", color = guide_legend(override.aes = list(linewidth = 3))) +
      scale_color_manual(
        values = color_values,
        labels = c("BRCA1/2m Cancer Negative",
                   cancer_type,
                   pos)
      ) +
      scale_alpha_manual(values = c("negative" = 0.25, "postiive" = 0.75)) +
      facet_wrap(.~diag, nrow=1) +
      theme +
      scale_x_continuous(limits = c(-900, 900), expand = c(0, 0))
    
    # Capture the two target plots:
    if (site == "BRCA" && list == "breast" && cancer_type == "breast")  UEC_breast_plot  <- fig
    if (site == "OV"   && list == "ovarian" && cancer_type == "ovarian") UEC_ovarian_plot <- fig
    
    ggsave(file.path(outdir, paste0("03_griffin_nucleosome_curves_TCGA_by_cancer_type_", site ,"_", cancer_type, ".pdf")), fig, width = 14, height = 7)
    plot_name <- paste0(site, "_plot")
    assign(plot_name, fig)
    
    # Print completion message for debugging
    cat("Completed processing for site:", site, "\n")
    
  } 
}

library(cowplot)
library(ggplot2)

# 1) Define the breast‐plot relabel map:
breast_labels <- c(
  "BRCA1/2m Cancer Negative"                     = "Cancer Negative",
  "BRCA1/2m Cancer Positive\n (Excluding breast)" = "Cancer Positive\n (Excluding Matched Cancer)",
  "breast"                                        = "Matched Cancer",
  "Healthy Control"                               = "Healthy Control"
)

# 2) Build UEC_breast_plot_fixed:
UEC_breast_plot_fixed <- UEC_breast_plot +
  scale_color_manual(values = c(
    "BRCA1/2m Cancer Negative"                     = "#76EE00",
    "BRCA1/2m Cancer Positive\n (Excluding breast)" = "#EE0000",
    "breast"                                        = "#FFA500",
    "Healthy Control"                               = "#2ca02c"
  )) +
  facet_wrap(
    ~ diag,
    nrow     = 1,
    labeller = labeller(diag = breast_labels)
  ) +
  labs(
    title = "Nucleosome Accesbility Coverage by Cancer Type",
    y     = "Breast Cancer Coverage"
  ) +
  theme(
    legend.position = "none",
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.line.x     = element_blank(),
    plot.title      = element_text(size = 16),
    strip.text      = element_text(size = 12),
    axis.title.y    = element_text(size = 14),
    axis.text.y     = element_text(size = 12)
  )

# 3) Define the ovarian‐plot relabel map:
ovarian_labels <- c(
  "BRCA1/2m Cancer Negative"                        = "Cancer Negative",
  "BRCA1/2m Cancer Positive\n (Excluding ovarian)"   = "Cancer Positive\n (Excluding Matched Cancer)",
  "ovarian"                                           = "Matched Cancer",
  "Healthy Control"                                   = "Healthy Control"
)

# 4) Make sure “diag” in UEC_ovarian_plot$data is factored in that same order:
UEC_ovarian_plot$data$diag <- factor(
  UEC_ovarian_plot$data$diag,
  levels = c(
    "BRCA1/2m Cancer Negative",
    "BRCA1/2m Cancer Positive\n (Excluding ovarian)",
    "ovarian",
    "Healthy Control"
  )
)

# 5) Build UEC_ovarian_plot_fixed:
UEC_ovarian_plot_fixed <- UEC_ovarian_plot +
  scale_color_manual(values = c(
    "BRCA1/2m Cancer Negative"                        = "#76EE00",
    "BRCA1/2m Cancer Positive\n (Excluding ovarian)"   = "#EE0000",
    "ovarian"                                           = "#FFA500",
    "Healthy Control"                                   = "#2ca02c"
  )) +
  facet_wrap(
    ~ diag,
    nrow     = 1,
    labeller = labeller(diag = ovarian_labels)
  ) +
  labs(
    y = "Ovarian Cancer Coverage"
  ) +
  theme(
    legend.position   = "none",
    plot.title        = element_blank(),
    strip.text        = element_blank(),
    strip.background  = element_blank(),
    axis.title.y      = element_text(size = 14),
    axis.text.y       = element_text(size = 14),
    axis.title.x      = element_text(size = 12),
    axis.text.x       = element_text(size = 10)
  )

# 6) Stack them in one column:
combined <- plot_grid(
  UEC_breast_plot_fixed,
  UEC_ovarian_plot_fixed,
  ncol        = 1,
  align       = "v",
  axis        = "lr",
  rel_heights = c(1, 1)
)

# 7) Adjust font sizes and save to PDF:
combined_final <- combined +
  theme(
    strip.text   = element_text(size = 14),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12)
  )

ggsave(
  filename = file.path(outdir, "03_griffin_breast_ovarian_combined.pdf"),
  plot     = combined_final,
  width    = 9.5,
  height   = 6,
  units    = "in"
)