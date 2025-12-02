library(tidyverse)
library(plyr)
library(GeneCycle)
library(matrixStats)
library(here)

process_griffin <- function(analysis_name, base_dir) {
  # Find all coverage files (recursive)
  all_files <- list.files(
    path = base_dir,
    pattern = "\\.coverage\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Exclude failed samples
  all_files <- all_files[!grepl("CHARMQ_0788_Pl_T_PG_T-788|CHARMQ_0207_Pl_T_PG_T-207", all_files)]
  
  # Output path
  outpath <- here::here("data", "griffin", "griffin_coverage")
  
  # Extract project (TCGA/custom) and site (BRCA, etc.)
  paths_df <- data.frame(path = all_files, stringsAsFactors = FALSE) %>%
    mutate(
      path_parts = strsplit(path, "/"),
      project = map_chr(path_parts, ~ .x[which(.x == "griffin") + 1]),
      site    = map_chr(path_parts, ~ .x[which(.x == "griffin") + 2])
    ) %>%
    dplyr::select(-path_parts) %>%
    distinct(project, site)
  
  # Loop through each project and site
  for (i in seq_len(nrow(paths_df))) {
    project <- paths_df$project[i]
    site <- paths_df$site[i]
    
  
    
    # Skip processing for excluded sites
    if (site %in% c("all_sites", "ZNF22", "TBX5")) {
      next
    }
    
    # Define the paths for this project and site
    current_paths <- all_files[
      grepl(paste0("/griffin/", project, "/", site, "/"), all_files)
    ]
    
    # Define directories and output folders dynamically
    path <- file.path(outpath, analysis_name, project, "coverage")
    outdir <- file.path(outpath, analysis_name, project)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # Create a loop to extract files and data from all samples for this project and site
    datalist <- lapply(current_paths, function(x) { read.delim(file = x, header = TRUE) })
    griffin <- do.call(rbind, datalist)
    
    # Keep columns for summary
    keep <- c(
      "GC_corrected_total_read_ends", "GC_corrected_total_read_starts",
      "background_normalization", "endpoint", "number_of_sites", "sample", 
      "site_name", "total_read_ends", "total_read_starts"
    )
    
    summary <- griffin[griffin$GC_correction == "GC_corrected", 
                       colnames(griffin) %in% keep]
    

    
    # Subset and process raw and corrected data
    raw <- griffin[griffin$GC_correction == "none", ]
    row.names(raw) <- raw$sample
    raw <- raw[, !(colnames(raw) %in% keep)]
    raw <- raw[, !(colnames(raw) %in% c("GC_correction"))]
    raw <- as.data.frame(t(raw))
    

    
    # Check header to extract distance
    if (file.exists(current_paths[1])) {
      header <- read.delim(current_paths[1], header = FALSE)
      header <- header[1, ]
      header <- header %>% dplyr::select(where(is.numeric))
      header <- t(header)
      raw$distance <- header
      raw <- raw[order(raw$distance), ]
    } else {
      stop("File does not exist: ", current_paths[1])
    }
    

    
    corrected <- griffin[griffin$GC_correction == "GC_corrected", ]
    row.names(corrected) <- corrected$sample
    corrected <- corrected[, !(colnames(corrected) %in% keep)]
    corrected <- corrected[, !(colnames(corrected) %in% c("GC_correction"))]
    corrected <- as.data.frame(t(corrected))
    corrected$distance <- header
    corrected <- corrected[order(corrected$distance), ]
    
    # Compute statistics
    stats <- corrected[, !(colnames(corrected) %in% "distance")]
    mean <- colMeans(stats)
    midpoint <- colMeans(stats[rownames(stats) %in% c("X.30", "X.15", "X0", "X15", "X30"), ])
    fft <- GeneCycle::periodogram(stats)[["spec"]]
    fft <- colMaxs(fft)
    stats <- rbind(mean, midpoint, fft)
    features <- c(paste0(site, "_coverage"), paste0(site, "_midpoint"), paste0(site, "_amp"))
    stats <- cbind(features, stats)
    
    # Write output
    write.table(summary, file.path(outdir, paste0(project, "_griffin_summary_", site, ".txt")),
                row.names = FALSE, sep = "\t")
    write.table(raw, file.path(outdir, paste0(project, "_griffin_raw_", site, ".txt")),
                row.names = FALSE, sep = "\t")
    write.table(corrected, file.path(outdir, paste0(project, "_griffin_corrected_", site, ".txt")),
                row.names = FALSE, sep = "\t")
    write.table(stats, file.path(outdir, paste0(project, "_griffin_features_", site, ".txt")),
                row.names = TRUE, sep = "\t")
    
    cat("Completed processing for site:", site, "and wrote to", outdir, "\n")
  }
  
}