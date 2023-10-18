# qPCR_setup_pooling
Code for manual qPCR setup and amplicon pooling for eDNA.

## Setup
Install tidyverse, readxl, and openxlsx. This can be done in R studio with: 
```
install.packages("tidyverse")
install.packages("readxl")
install.packages("openxlsx")
```

### This doesn't work yet because it's a private repo
```
library(devtools)
install_github("MarcieAyad/qPCR_setup_pooling")
```

## meta_to_plates() function
### Running meta_to_plates()
This is the minimum arguments required to run the meta_to_plates() function:
```
input_file  <- "path/to/input.xlsx"
output_file <- "path/to/output.xlsx"
assays      <- c("16S", "MiFish")
run         <- "run_1"

meta_to_plates(
    input_file,
    output_file,
    assays,
    run
)
```

This is all the arguments that can be used when running the meta_to_plates() function:
```
input_file      <- "path/to/input.xlsx"
output_file     <- "path/to/output.xlsx"
assays          <- c("16S", "MiFish")
run             <- "run_1"
plate_width     <- 12
plate_height    <- 8
controls        <- c("ITC", "NTC")
control_pattern <- "WC|DI|EB|BC|NTC|ITC|Cont"

meta_to_plates(
    input_file,
    output_file,
    assays,
    run,
    plate_width,
    plate_height,
    controls,
    control_pattern
)
```

### meta_to_plates() arguments
- `input_file`: A path to an Excel file in the format explained below in the 'meta_to_plates() input' section.
- `output_file`: The path for your output Excel file.
- `assays`: A vector of assays. These should match the assay names found in the input file.
- `run`: The name of the run. This should match the sequencing run name found in the input file. 
- `plate_width`: The number of columns per plate. Default = 12.
- `plate_height`: The number of rows per plate. Default = 8. Max allowed = 13. 
- `controls`: A vector of control samples that will be added at the end of each plate. Default = c("ITC", "NTC").
- `control_pattern`: A string that will be used to flag samples as control samples. Default = "WC|DI|EB|BC|NTC|ITC|Cont". This default value means that any sample containing "WC", "DI", "EB", "BC", "NTC", "ITC", or "Cont" will be flagged as a control sample.

### meta_to_plates() input
An example input file can be viewed at `test_data/fake_meta.xlsx`

The input Excel file should contain a metadata sheet and one index sheet for each assay.
- `metadata`: Should have the columns `sample_id` and `sequencing_run`.
  - `sample_id`: This column will be used to name your samples and to flag samples as control samples.
  - `sequencing_run`: This column will be used to use only samples with the same run as the `run argument`.
- `${assay}_index`: Replace `${assay}` with the name of your assay. Should have the columns `primer_num`, `primer_seq`, `tags`, and `fw_rv`.
  - `primer_num`: Each value here should be unique. These values will be used as column/row names for your plates.
  - `primer_seq`: The primer sequence. This will be added to the output file.
  - `tags`: The sequence used to for demultiplexing samples.
  - `fw_rv`: Use this column to indicate if the primer is a fw or rv primer. Values can be `fw` or `rv`.

### meta_to_plates() output
An example output file can be viewed at `test_data/results.xlsx`

The output Excel file will have four sheets; `plates`, `big_plates`, `metadata`, and `position_df`.
- `plates`: Will be all the plates with their.
- `big_plates`: Will be bigger versions of the plates with 3 replicates and a pool replicate added to each sample.
- `metadata`: Your samples will be duplicated for each assay. The metadata will now contain demultiplex, plate, and well information.
- `position_df`: This will have your samples duplicated for each replicate in the big plates. This sheet is needed for the `plates_to_pooling()` function.

## plates_to_pooling() function
Not finished yet

### Running plates_to_pooling()
### plates_to_pooling() arguments
### plates_to_pooling() input
### plates_to_pooling() output
