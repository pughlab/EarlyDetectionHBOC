# --- 1. CRAN packages ---
cran_pkgs <- c(
  "glue", "gbm", "broom", "ggsignif", "car", "carData",
  "combinat", "ggforce", "caret", "lattice", "lemon",
  "signal", "ggformula", "ggridges", "scales", "pracma",
  "gsignal", "GeneCycle", "fdrtool", "longitudinal",
  "corpcor", "MASS", "pROC", "psych", "RColorBrewer",
  "ggh4x", "reshape2", "gridExtra", "ggpubr", "cowplot",
  "circlize", "matrixStats", "data.table", "plyr",
  "lubridate", "forcats", "purrr", "readr", "tidyr",
  "tibble", "ggplot2", "Polychrome", "readxl",
  "stringr", "dplyr", "here", "rmarkdown",
  # dependencies already seen in sessionInfo
  "magick", "survival", "parallelly", "sass", "bslib",
  "future", "clue", "digest", "colorspace", "Matrix",
  "rprojroot", "textshaping", "labeling", "randomForest",
  "timechange", "polyclip", "abind", "proxy", "withr",
  "lava", "rjson", "scatterplot3d", "ModelMetrics",
  "nlme", "cluster", "generics", "recipes", "gtable",
  "labelled", "tzdb", "class", "hms", "foreach",
  "pillar", "splines", "tweenr", "xfun", "hardhat",
  "mosaicCore", "timeDate", "stringi", "yaml",
  "evaluate", "codetools", "kernlab", "cli", "rpart",
  "systemfonts", "munsell", "Rcpp", "globals", "png",
  "gower", "listenv", "ipred", "e1071", "prodlim",
  "crayon", "GetoptLong", "rlang", "mnormt"
)

# --- 2. Bioconductor packages --
bioc_pkgs <- c(
  "ComplexHeatmap",
  "S4Vectors",
  "IRanges",
  "BiocGenerics"
)

# --- 3. Install CRAN packages if missing ---
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# --- 4. Install Bioconductor packages if missing ---
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_bioc) > 0) {
  if (!"BiocManager" %in% installed.packages()[, "Package"]) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc, ask = FALSE)
}

# --- 5. Load everything ---
all_pkgs <- c(cran_pkgs, bioc_pkgs)

invisible(lapply(all_pkgs, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))

message("All CRAN and Bioconductor packages successfully loaded.")