
# qPCR_setup_pooling

Code for manual qPCR setup and amplicon pooling for eDNA.

## Setup

Install dependencies. This can be done in R studio with:

```{R}
install.packages("tidyverse")
install.packages("readxl")
install.packages("openxlsx")
install.packages("plyr")
install.packages("gtools")
install.packages("ggplate")
```

## Script templates

There is a template file in the R folder called `run_meta_to_plates.R` that can be used to run the meta_to_plates() function. The full script to run the plates-to-pooling functionality can be found in `run_plates_to_pooling.R`. The script to merge qPCR QC metrics can be found in `run_concat_qPCR_data.R`.

## run_concat_qPCR_data.R

### Running run_concat_qPCR_data

These are the variables you may need to change in the `run_concat_qPCR_data.R` script:

```{R}
input_dir  <- "test_data/input/cp_epf_tm_files/"
output_dir <- "test_data/output/"
assays     <- c("16S", "MiFish")
```

### run_concat_qPCR_data variables

- `input_dir`: Directory with qPCR QC data. Filenames should be in a specific format mentioned below.
- `output_dir`: Directory where you would like the output file. The output file will contain the Cp, EPF, and Tm data in one file.
- `assays`: A vector of assays. These should match the assay names found in the input file names.

### run_concat_qPCR_data input

- `input_dir`: An example of the qpcr directory can be viewed at `test_data/cp_epf_tm_files`
  - Each file in the qPCR directory should have names that end in `_$assay_$plate_Cp.txt`, `_$assay_$plate_EPF.txt`, or `_$assay_$plate_Tm.txt`. An example name would be `20240521_SWWA_V7_16S_Plate1_SH_Cp.txt`.
  - These txt files should be tab seperated files matching the format found in those test files.

## Running meta_to_plates()

This is the minimum arguments required to run the meta_to_plates() function:

```{R}
source("R/meta_to_plates.R")
metadata    <- "path/to/metadata.xlsx"
index_file  <- "path/to/indexes.xlsx"
output_file <- "path/to/output.xlsx"
assays      <- c("16S", "MiFish")

meta_to_plates(
    input_file,
    output_file,
    assays
)
```

Note: This assumes the R folder with all the scripts is in your current working directory.

This is all the arguments that can be used when running the meta_to_plates() function:

```{R}
metadata        <- "path/to/metadata.xlsx"
index_file      <- "path/to/indexes.xlsx"
output_file     <- "path/to/output.xlsx"
assays          <- c("16S", "MiFish")
run             <- "run_1"
plate_width     <- 12
plate_height    <- 8
controls        <- c("NTC", "ITC")
control_pattern <- "WC|DI|EB|BC|NTC|ITC|Cont|BL"
skip_plates     <- list(
  "16S" = 0,
  "MiFish" = 0
)
skip_samples    <- list(
  "16S" = 0,
  "MiFish" = 0
)
strategy        <- "UC"

meta_to_plates(
    input_file,
    output_file,
    assays,
    run,
    plate_width,
    plate_height,
    controls,
    control_pattern,
    skip_plates,
    skip_samples,
    strategy
)
```

### meta_to_plates() arguments

- `metadata`: A path to an Excel file that contains your metadata. Make sure you have a 'samp_name' column with the unique identifiers for your samples. This is compatible with FAIRe formatted metadata. This also accepts multiple Excel files if you provide a vector of paths.
- `input_file`: A path to an Excel file in the format explained below in the 'meta_to_plates() input' section.
- `output_file`: The path for your output Excel file.
- `assays`: A vector of assays. These should match the assay names found in the index file.
- `run`: The name of the run. This should match the 'sequencing_run' name found in the metadata. Samples will only be used if they match the 'run'. This feature won't be used if the metadata doesn't have a 'sequencing_run' column.
- `plate_width`: The number of columns per plate. Default = 12.
- `plate_height`: The number of rows per plate. Default = 8. Max allowed = 13.
- `controls`: A vector of control samples that will be added at the end of each plate. Default = c("NTC", "ITC").
- `control_pattern`: A string that will be used to flag samples as control samples. Default = "WC|DI|EB|BC|NTC|ITC|Cont|BL". This default value means that any sample containing "WC", "DI", "EB", "BC", "NTC", "ITC", "Cont", or "BL" will be flagged as a control sample.
- `skip_plates`: How many plates would you like to skip? Allows different values for each assay.
- `skip_samples`: How many samples would you like to skip? Allows different values for each assay. Can be used with skip_plates (e.g., you can skip 2 plates and 10 samples).
- `strategy`: Can be "UDI" for unique dual-index, or "UC" for unique combinatorial

### meta_to_plates() input

An example input file can be viewed at `test_data/AB_V12_V9_metadata.xlsx`
There is also a template you can follow at `test_data/R_input_Template.xlsx`

The index Excel file should contain one index sheet for each assay.

- `${assay}_index`: Replace `${assay}` with the name of your assay. Should have the columns `primer_num`, `tags`, and `fw_rv`.
  - `primer_num`: Each value here should be unique. These values will be used as column/row names for your plates.
  - `tags`: The sequence used to for demultiplexing samples.
  - `fw_rv`: Use this column to indicate if the primer is a fw or rv primer. Values can be `fw` or `rv`.

### meta_to_plates() output

An example output file can be viewed at `test_data/output.xlsx`

The output Excel file will have four sheets; `plates`, `big_plates`, `metadata`, and `position_df`.

- `plates`: Will be all the plates with their.
- `big_plates`: Will be bigger versions of the plates with 3 replicates and a pool replicate added to each sample.
- `metadata`: Your samples will be duplicated for each assay. The metadata will now contain demultiplex, plate, and well information.
- `position_df`: This will have your samples duplicated for each replicate in the big plates. This sheet is needed for the `run_plates_to_pooling()` script.
- `samplesheet_${assay}`

## plates_to_pooling() function

### Running run_plates_to_pooling

These are the variables you may need to change in the `run_plates_to_pooling.R` script:

```{R}
input_file   <- "test_data/output/output_df.xlsx"
qpcr_dir     <- "test_data/input/qPCR_test_data"
output_dir   <- "test_data/output/"
plate_width  <- 12
plate_height <- 8
assay        <- c("16S", "MiFish")
suffix       <- ""
```

### run_plates_to_pooling variables

- `input_file`: A path to the Excel file that was the output of the meta_to_plates() function.
- `qpcr_dir`: Directory with qPCR data. Filenames should be in a specific format mentioned below.
- `output_dir`: Directory where you would like all output files. The output files are explained in more detail below.
- `plate_width`: The number of columns per plate. Default = 12. This should be the same value you used in meta_to_plates().
- `plate_height`: The number of rows per plate. Default = 8. Max allowed = 13. This should be the same value you used in meta_to_plates().
- `assay`: The assays that you're working with.
- `suffix`: What suffix would you like for your file names.

### run_plates_to_pooling input

- `input_file`: This is the output file from the meta_to_plates() function. An example can be viewed at `test_data/results.xlsx`
- `qpcr_dir`: An example of the qpcr directory can be viewed at `test_data/qPCR_test_data`
  - Each file in the qPCR directory should have names that end in `_$assay_$plate.txt`. An example name would be `20230901_Extraction_16S_Plate1.txt`.
  - These txt files should be tab separated files with the columns `Experiment_name`, `Position`, `Sample`, `EPF`, `Cp`, `Tm1`, `Tm2`.

### run_plates_to_pooling discarding samples

During the running of the run_plates_to_pooling script, there is a section with the header `Output rxns_to_check.csv`.
This section will create a file in your `output_dir` called `rxns_to_check.csv`.
You can manually change the `discard` column in this file to choose replicates that you would like to keep or discard.
The values in the `discard` column should be either `KEEP` or `DISCARD`.
This file will be imported back into the script and your changes to the `discard` column will be merged.

### run_plates_to_pooling output

The output directory will contain multiple .csv and .pdf files.

- `rxns_to_check.csv`: This file is explained in more detail above.

- `reps_to_discard.csv`: This file has information on replicates that have been flaged for discarding.
  - `assay`: This column will have the assay of the discard replicate.
  - `plate_number`: The plate number of the discard replicate.
  - `pos`: The position on the plate of the discard replicate.
  - `sample_replicate`: The sample name and the replicate number of the discard replicate.
  - `discard`: The flag DISCARD.
  - `sample`: The sample name.

- `biomek_pooling_workbook_$n.csv`: There will be one of these files for every four plates.
  - `assay`: The assay.
  - `plate_number`: The plate number.
  - `SourcePosition`: The source position starting from P11 and going down to P8.
  - `SourceWell`: The source well.
  - `Volume`: The volume in uL.
  - `DestinationPosition`: The destination position.
  - `DestinationWell`: The destination well.

- `raw_$assay_$plate.pdf`: Plate plots before cleaning.

- `clean_$assay_$plate.pdf`: Plate plots after cleaning.

- `$assay_$plate_$sampletype_minipool_plot.pdf`: Sample type specific minipool plots.
