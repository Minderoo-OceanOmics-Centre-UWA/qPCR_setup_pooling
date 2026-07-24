source("R/meta_to_plates.R")
metadata    <- "test_data/input/AB_V12_V9_metadata.xlsx" # Use a vector of paths of you have multiple metadata files
index_file  <- "input/index_Template.xlsx"
output_dir  <- "output/"
assays      <- c("16SFishD", "MarVer1", "MiFishUE2")
prefix      <- "" 

plates_to_skip <- list(
  "16SFishD" = 0,
  "MarVer1" = 0,
  "MiFishUE2" = 0
)
samples_to_skip <- list(
  "16SFishD" = 0,
  "MarVer1" = 0,
  "MiFishUE2" = 0
)

# strategy can be 'UDI' for unique dual-index, or 'UC' for unique combinatorial

meta_to_plates(
  metadata,
  index_file,
  output_file,
  assays,
  skip_plates = plates_to_skip,
  skip_samples = samples_to_skip,
  strategy = "UC",
  prefix = prefix
)

