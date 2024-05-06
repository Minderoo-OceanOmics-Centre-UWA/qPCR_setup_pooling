source("R/meta_to_plates.R")
input_file  <- "test_data/input/AB_V12_V9_metadata.xlsx"
output_file <- "test_data/output/output_df.xlsx"
assays      <- c("16S", "MiFish")
run         <- "run1"

meta_to_plates(
  input_file,
  output_file,
  assays,
  run
)
