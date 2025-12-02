library(here)
HBOC <- read.table(here("raw_data", "griffin", "2025-10-24_HBOC_WGS_downsampled_bwa_bams_16102025_nucleosome_accessibility_distances.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

healthy <- read.table(here("raw_data", "griffin", "2025-10-30_frag_pipeline_HBC_27102025_nucleosome_accessibility_distances.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Save the combined healthy dataset
saveRDS(healthy, file = here("raw_data", "griffin", "nucleosome_accessibility_distances_HBC.rds"))

# Save the HBOC dataset separately
saveRDS(HBOC, file = here("raw_data", "griffin", "nucleosome_accessibility_distances_HBOC.rds"))