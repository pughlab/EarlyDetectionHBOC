#!/usr/bin/env Rscript

# Set global seed for reproducibility
set.seed(42)

library(rmarkdown)
library(here)

cat("\n========================================\n")
cat("       HBOC Figure Generation Pipeline\n")
cat("========================================\n\n")


### ================ Set up ================

download_raw_data <- function(url, dest_dir = "raw_data") {
  # Temporary path for the zip
  zip_path <- "raw_data.zip"
  
  message("Downloading raw_data.zip from Zenodo...")
  utils::download.file(url, zip_path, mode = "wb")
  
  message("Unzipping into project root...")
  utils::unzip(zip_path, exdir = ".")
  
  # Remove zip and __MACOSX if created
  unlink(zip_path)
  unlink("__MACOSX", recursive = TRUE)
  
  # If the archive contains raw_data/, we may want to move it
  if (dir.exists("raw_data/raw_data")) {
    file.rename("raw_data/raw_data", "raw_data_new")
    unlink("raw_data", recursive = TRUE)
    file.rename("raw_data_new", "raw_data")
  }
  
  message("Raw data download complete.")
}

zenodo_url <- "https://zenodo.org/" 

if (!dir.exists("raw_data")) {
  options(timeout = max(1800, getOption("timeout")))
  download_raw_data(zenodo_url)
}

# List of directories you want to ensure exist
dirs <- c(
  "data",
  "figures",
  "HTML",
  # add subfolders if you want them auto-created too
  "data/cfMeDIP",
  "data/classifier",
  "data/cohort",
  "data/coverage",
  "data/DNA_extraction",
  "data/fragment_length",
  "data/fragment_ratio",
  "data/griffin",
  "data/HBOC_pipeline_output",
  "data/ichorCNA",
  "data/longitudinal",
  "data/metrics",
  "data/mutations",
  "data/nucleosome_peaks",
  "data/statistics",
  "figures/cfMeDIP",
  "figures/classifier",
  "figures/cohort",
  "figures/fragment_length",
  "figures/fragment_ratio",
  "figures/griffin",
  "figures/ichorCNA",
  "figures/longitudinal",
  "figures/metrics",
  "figures/mutations",
  "figures/nucleosome_peaks",
  "figures/Statistics"
)

# Create dir if they don’t exist
for (d in dirs) {
  full_path <- here::here(d)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE)
    message("Created: ", full_path)
  }
}

# Create HTML output directory if it doesn't exist
html_dir <- here::here("HTML")
if (!dir.exists(html_dir)) dir.create(html_dir, recursive = TRUE)

# =============================================================================
# Mutations Oncoprint
# =============================================================================
cat("\n========== Mutations Oncoprint ==========\n")

render(
  input       = here::here("scripts", "mutations", "02_figure_mutations_oncoplot_to_publish.Rmd"),
  output_file = here::here("HTML", "02_figure_mutations_oncoplot_to_publish.html"),
  clean       = TRUE,
  params      = list(ichorCNA_estimate_table = "HBOC_ichorCNA_estimates.tsv", HBC_ichorCNA_estimate_table = "HBC_ichorCNA_estimates.tsv")
)

cat("\nMutation oncoprint generated.\n")

# =============================================================================
# ichorCNA Segmentation Processing
# =============================================================================
cat("\n========== ichorCNA Segmentation ==========\n")

render(
  input       = here::here("scripts", "ichorCNA", "01_processing_ichorCNA_segs.Rmd"),
  output_file = here::here("HTML", "01_processing_ichorCNA_segs.html"),
  clean       = TRUE, 
  params      = list(raw_data_input_dir = "HBOC_ichorCNA")
)

cat("\nichorCNA segmentation processing completed.\n")

# =======================================================
# Figure 1 ichorCNA
# =======================================================
cat("\n========== Figure 1 ichorCNA ==========\n")

render(
  input       = here::here("scripts", "ichorCNA", "02_figure_ichorCNA_correctedDepth_to_publish.Rmd"),
  output_file = here::here("HTML", "02_figure_ichorCNA_correctedDepth_to_publish.html"),
  clean       = TRUE,
  params      = list(ichorCNA_estimate_table = "HBOC_ichorCNA_estimates.tsv")
)

cat("\nFigure 1 ichorCNA generated.\n")


# =============================================================================
# Fragmentomics Processing (HBOC + Control)
# =============================================================================
cat("\n========== Fragmentomics Processing ==========\n")

# =============================
# Collect Fragmentomics - HBOC
# =============================
render(
  input       = here::here("scripts", "HBOC_pipeline_output", "01_HBOC_pipeline_all_output.Rmd"),
  output_file = here::here("HTML", "01_HBOC_pipeline_all_output_HBOC.html"),
  params      = list(which_project = "HBOC"),
  clean       = TRUE
)

# =============================
# Collect Fragmentomics - Control
# =============================
render(
  input       = here::here("scripts", "HBOC_pipeline_output", "01_HBOC_pipeline_all_output.Rmd"),
  output_file = here::here("HTML", "01_HBOC_pipeline_all_output_control.html"),
  params      = list(which_project = "control"),
  clean       = TRUE
)

# =========================================================================
# Fragment Length - 600bp
# =========================================================================
cat("\n========== Fragment Length 600bp ==========\n")

render(
  input       = here::here("scripts", "fragment_length", "01_figure_fragment_length_600bp.Rmd"),
  output_file = here::here("HTML", "01_figure_fragment_length_600bp.html"),
  clean       = TRUE
)
cat("\nFragment length 600bp figure generated.\n")

# =============================
# Fragment Length - Mono/Di/Tri Nucleosomes
# =============================
cat("\n========== Fragment Length Mono/Di/Tri Nucleosomes ==========\n")

render(input = here::here("scripts","fragment_length","02_figure_fragment_length_mono_di_tri_nucl_pos_vs_neg.Rmd"),
       output_file = here::here("HTML","02_figure_fragment_length_mono_di_tri_nucl_pos_vs_neg.html"),
       clean = TRUE)

cat("\nFragment Length Mono/Di/Tri Nucleosomes figure generated.\n")

# =============================
# Fragment Length - Mono/Di/Tri Nucleosome
# =============================
cat("\n========== Fragment Length Mono/Di/Tri Nucleosome ==========\n")
render(input = here::here("scripts","fragment_length","03_figure_fragment_length_mono_di_tri_nucl.Rmd"),
       output_file = here::here("HTML","03_figure_fragment_length_mono_di_tri_nucl.html"),
       clean = TRUE)

# =============================
# Fragment Ratio
# =============================
cat("\n========== Fragment Ratio ==========\n")

# 01 - Processing Fragment Ratio Statistics
render(
  input       = here::here("scripts", "fragment_ratio", "01_processing_fragment_ratio_statistics.Rmd"),
  output_file = here::here("HTML", "01_processing_fragment_ratio_statistics.html"),
  clean       = TRUE
)

# 02 - Heatmap to Publish
render(
  input       = here::here("scripts", "fragment_ratio", "02_figure_fragment_ratio_heatmap_to_publish.Rmd"),
  output_file = here::here("HTML", "02_figure_fragment_ratio_heatmap_to_publish.html"),
  clean       = TRUE,
  params      = list(ichorCNA_estimate_table = "HBOC_ichorCNA_estimates.tsv", HBC_ichorCNA_estimate_table = "HBC_ichorCNA_estimates.tsv")
)

# 03 - Profiles with Box Plot to Publish
render(
  input       = here::here("scripts", "fragment_ratio", "03_figure_fragment_ratio_profiles_with_box_plot_to_publish.Rmd"),
  output_file = here::here("HTML", "03_figure_fragment_ratio_profiles_with_box_plot_to_publish.html"),
  clean       = TRUE
)

# =============================
# Nucleosome Peaks 
# =============================

cat("\n========== Nucleosome Peaks by Source ==========\n")

render(
  input       = here::here("scripts", "nucleosome_peaks", "01_figure_nucleosome_peaks_by_source_to_publish.Rmd"),
  output_file = here::here("HTML", "01_figure_nucleosome_peaks_by_source_to_publish.html"),
  clean       = TRUE
)

cat("\n========== Nucleosome Peaks Score Calculation ==========\n")

render(
  input       = here::here("scripts", "nucleosome_peaks", "02_processing_nucleosome_peaks_healthy_correlation.Rmd"),
  output_file = here::here("HTML", "02_processing_nucleosome_peaks_healthy_correlation.html"),
  clean       = TRUE
)

# ======================================================
# ===================  GRIFFIN ANALYSIS  ===================
# ==========  Full Pipeline: HBOC + Control  ==========
# ======================================================

cat("Starting Griffin processing pipeline...\n")

# Step 0: Load Functions
source(here::here("scripts", "griffin", "01_processing_griffin.R"))

# Step 1: Process HBOC samples
cat("Step 1: Processing HBOC samples...\n")
process_griffin(
  analysis_name = "01_processing_griffin_HBOC",
  base_dir = here::here(
    "raw_data", "HBOC_pipeline_output", 
    "HBOC_organized_outputs", "griffin"
  )
)

# Step 2: Process Control samples
cat("Step 2: Processing Control samples...\n")
process_griffin(
  analysis_name = "01_processing_griffin_control",
  base_dir = here::here(
    "raw_data", "HBOC_pipeline_output", 
    "HBC_organized_outputs", "griffin"
  )
)

# Step 3: Generate Coverage Figures
cat("Step 3: Generating coverage figures...\n")
source(here::here("scripts", "griffin", "02_figure_griffin_coverage_by_source.R"))
source(here::here("scripts", "griffin", "03_figure_griffin_coverage_by_cancer_to_publish.R"))

# Step 4: Midpoint Amplitude Dot Plot (RMarkdown)
cat("Step 4: Rendering midpoint amplitude dot plot...\n")
rmarkdown::render(
  input       = here::here("scripts", "griffin", "04_figure_griffin_midpoint_amplitude_dot_plot.Rmd"),
  output_file = here::here("HTML",   "04_figure_griffin_midpoint_amplitude_dot_plot.html"),
  clean       = TRUE
)

# Step 5: Griffin Score Calculation
cat("Step 5: Calculating Griffin scores...\n")
source(here::here("scripts", "griffin", "05_processing_statistics_griffin_healthy_correlation.R"))

# Step 6: Griffin Scoring Heatmap (RMarkdown)
cat("Step 6: Rendering scoring heatmap...\n")
rmarkdown::render(
  input       = here::here("scripts", "griffin", "06_figure_griffin_scoring_heatmap_matched_cancer_to_publish.Rmd"),
  output_file = here::here("HTML",   "06_figure_griffin_scoring_heatmap_matched_cancer_to_publish.html"),
  clean       = TRUE
)

cat("Griffin processing pipeline completed successfully!\n")

# ==========================
# Statistics
# ==========================

# 01_processing_statistics.Rmd
rmarkdown::render(
  input       = here::here("scripts", "statistics", "01_processing_statistics.Rmd"),
  output_file = here::here("HTML", "01_processing_statistics.html"),
  clean       = TRUE,
  params      = list(ichorCNA_estimate_table = "HBOC_ichorCNA_estimates.tsv", HBC_ichorCNA_estimate_table = "HBC_ichorCNA_estimates.tsv")
)

# 02_processing_statistics_integrated_cancer_scores.Rmd
rmarkdown::render(
  input       = here::here("scripts", "statistics", "02_processing_statistics_integrated_cancer_scores.Rmd"),
  output_file = here::here("HTML", "02_processing_statistics_integrated_cancer_scores.html"),
  clean       = TRUE
)

# 03_figure_statistics_integrated_barplots_by_cancer_status_to_publish.Rmd
rmarkdown::render(
  input       = here::here("scripts", "statistics", "03_figure_statistics_integrated_barplots_by_cancer_status_to_publish.Rmd"),
  output_file = here::here("HTML", "03_figure_statistics_integrated_barplots_by_cancer_status_to_publish.html"),
  clean       = TRUE
)

# 04_figure_statistics_integrated_barplots_by_LIB_ID_to_publish.Rmd
rmarkdown::render(
  input       = here::here("scripts", "statistics", "04_figure_statistics_integrated_barplots_by_LIB_ID_to_publish.Rmd"),
  output_file = here::here("HTML", "04_figure_statistics_integrated_barplots_by_LIB_ID_to_publish.html"),
  clean       = TRUE
)

# ==========================
# Longitudinal
# ==========================

# 01_processing_cohort_longitudinal_info.Rmd
rmarkdown::render(
  input       = here::here("scripts", "longitudinal", "01_processing_cohort_longitudinal_info.Rmd"),
  output_file = here::here("HTML", "01_processing_cohort_longitudinal_info.html"),
  clean       = TRUE
)

# 02_figure_statistics_longitudinal_scores_4_cases.Rmd
rmarkdown::render(
  input       = here::here("scripts", "longitudinal", "02_figure_statistics_longitudinal_scores_4_cases.Rmd"),
  output_file = here::here("HTML", "02_figure_statistics_longitudinal_scores_4_cases.html"),
  clean       = TRUE
)

# 03_longitudinal_general_statistics.Rmd
rmarkdown::render(
  input       = here::here("scripts", "longitudinal", "03_longitudinal_general_statistics.Rmd"),
  output_file = here::here("HTML", "03_longitudinal_general_statistics.html"),
  clean       = TRUE
)

# ==========================
# Metrics
# ==========================

# 01_metrics_integrated_processing_export_data.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "01_metrics_integrated_processing_export_data.Rmd"),
  output_file = here::here("HTML", "01_metrics_integrated_processing_export_data.html"),
  clean       = TRUE
)

# 02_figure_metrics_confusion.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "02_figure_metrics_confusion.Rmd"),
  output_file = here::here("HTML", "02_figure_metrics_confusion.html"),
  clean       = TRUE
)

# 03_figure_metrics_score_significance_test_to_publish.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "03_figure_metrics_score_significance_test_to_publish.Rmd"),
  output_file = here::here("HTML", "03_figure_metrics_score_significance_test_to_publish.html"),
  clean       = TRUE
)

# 04_metrics_integration_score_info_for_cohort_to_publish.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "04_metrics_integration_score_info_for_cohort_to_publish.Rmd"),
  output_file = here::here("HTML", "04_metrics_integration_score_info_for_cohort_to_publish.html"),
  clean       = TRUE
)

# 05_metrics_variants_copy_number.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "05_metrics_variants_copy_number.Rmd"),
  output_file = here::here("HTML", "05_metrics_variants_copy_number.html"),
  clean       = TRUE
)

# 06_metrics_individual_fragmentomics.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "06_metrics_individual_fragmentomics.Rmd"),
  output_file = here::here("HTML", "06_metrics_individual_fragmentomics.html"),
  clean       = TRUE
)

# 07_metrics_integrated_analysis.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "07_metrics_integrated_analysis.Rmd"),
  output_file = here::here("HTML", "07_metrics_integrated_analysis.html"),
  clean       = TRUE
)

# 08_metrics_processing_coverage_information.Rmd
rmarkdown::render(
  input       = here::here("scripts", "metrics", "08_metrics_processing_coverage_information.Rmd"),
  output_file = here::here("HTML", "08_metrics_processing_coverage_information.html"),
  clean       = TRUE
)

# ==========================
# Classifier
# ==========================

rmarkdown::render(
  input       = here::here("scripts", "classifier", "01_processing_classifier.Rmd"),
  output_file = here::here("HTML", "01_processing_classifier.html"),
  clean       = TRUE
)

# ==========================
# Cohort
# ==========================

rmarkdown::render(
  input       = here::here("scripts", "cohort", "01_processing_cohort_breakdown_with_cfMeDIP.Rmd"),
  output_file = here::here("HTML", "01_processing_cohort_breakdown_with_cfMeDIP.html"),
  clean       = TRUE
)

rmarkdown::render(
  input       = here::here("scripts", "cohort", "02_figure_cohort_swimmer_with_cfMeDIP.Rmd"),
  output_file = here::here("HTML", "02_figure_cohort_swimmer_with_cfMeDIP.html"),
  clean       = TRUE
)

rmarkdown::render(
  input       = here::here("scripts", "cohort", "03_figure_cohort_upset_plot_with_cfMeDIP.Rmd"),
  output_file = here::here("HTML", "03_figure_cohort_upset_plot_with_cfMeDIP.html"),
  clean       = TRUE
)

# ==========================
# cfMeDIP
# ==========================

# 01_cfMeDIP_combined_genomic_and_methylation.Rmd
rmarkdown::render(
  input       = here::here("scripts", "cfMeDIP", "01_cfMeDIP_combined_genomic_and_methylation.Rmd"),
  output_file = here::here("HTML", "01_cfMeDIP_combined_genomic_and_methylation.html"),
  clean       = TRUE
)

# 02_cfMeDIP_statistics.Rmd
rmarkdown::render(
  input       = here::here("scripts", "cfMeDIP", "02_cfMeDIP_statistics.Rmd"),
  output_file = here::here("HTML", "02_cfMeDIP_statistics.html"),
  clean       = TRUE
)

# 03_figure_cfMeDIP_barplots_by_cancer_status_w_methylation.Rmd
rmarkdown::render(
  input       = here::here("scripts", "cfMeDIP", "03_figure_cfMeDIP_barplots_by_cancer_status_w_methylation.Rmd"),
  output_file = here::here("HTML", "03_figure_cfMeDIP_barplots_by_cancer_status_w_methylation.html"),
  clean       = TRUE
)

# ==========================
# Coverage
# ==========================

# 01_coverage_information_WGS_TS_cfMeDIP.Rmd
rmarkdown::render(
  input       = here::here("scripts", "coverage", "01_coverage_information_WGS_TS_cfMeDIP.Rmd"),
  output_file = here::here("HTML", "01_coverage_information_WGS_TS_cfMeDIP.html"),
  clean       = TRUE
)

# ==========================
# DNA_extraction
# ==========================

# 01_metrics_DNA_extraction_HBOC_samples.Rmd
rmarkdown::render(
  input       = here::here("scripts", "DNA_extraction", "01_metrics_DNA_extraction_HBOC_samples.Rmd"),
  output_file = here::here("HTML", "01_metrics_DNA_extraction_HBOC_samples.html"),
  clean       = TRUE
)

# 02_metrics_DNA_extraction_control_samples.Rmd
rmarkdown::render(
  input       = here::here("scripts", "DNA_extraction", "02_metrics_DNA_extraction_control_samples.Rmd"),
  output_file = here::here("HTML", "02_metrics_DNA_extraction_control_samples.html"),
  clean       = TRUE
)

# =======================================================
# Completion
# =======================================================
cat("\n========================================\n")
cat("           All Figures Generated\n")
cat("========================================\n\n")