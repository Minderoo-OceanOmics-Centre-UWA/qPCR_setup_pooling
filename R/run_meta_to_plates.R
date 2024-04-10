source("R/meta_to_plates.R")
input_file  <- "test_data/fake_meta.xlsx"
output_file <- "test_data/test_output.xlsx"
assays      <- c("16S", "MiFish")
run         <- "run1"

meta_to_plates(
  input_file,
  output_file,
  assays,
  run
)
