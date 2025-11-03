source("R/meta_to_plates.R")
metadata    <- "test_data/input/AB_V12_V9_metadata.xlsx" # Use a vector of paths of you have multiple metadata files
index_file  <- "test_data/input/R_input_Template.xlsx"
output_file <- "test_data/output/output_df.xlsx"
assays      <- c("16S", "MiFishU", "MiFishE2", "COILeray")

plates_to_skip <- list(
  "16S" = 0,
  "MiFishU" = 0,
  "MiFishE2" = 0,
  "COILeray" = 0
)
samples_to_skip <- list(
  "16S" = 0,
  "MiFishU" = 0,
  "MiFishE2" = 0,
  "COILeray" = 0
)

# strategy can be 'UDI' for unique dual-index, or 'UC' for unique combinatorial

meta_to_plates(
  metadata,
  index_file,
  output_file,
  assays,
  skip_plates = plates_to_skip,
  skip_samples = samples_to_skip,
  strategy = "UC"
)

