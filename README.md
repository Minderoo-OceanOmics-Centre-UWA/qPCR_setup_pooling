# qPCR_setup_pooling
Code for manual qPCR setup and amplicon pooling for eDNA.

## Setup
Install tidyverse, readxl, and openxlsx. This can be done in R studio with: 
```
install.packages("tidyverse")
install.packages("readxl")
install.packages("openxlsx")
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
### Running plates_to_pooling()
This is the minimum arguments required to run the plates_to_pooling() function:
```
input_file <- "path/to/results.xlsx"
qpcr_dir   <- "path/to/qPCR_data"
output_dir <- "path/to/output_dir"

plates_to_pooling(
    input_file,
    qpcr_dir,
    output_dir
)
```

This is all the arguments that can be used when running the plates_to_pooing() function:
```
input_file   <- "path/to/results.xlsx"
qpcr_dir     <- "path/to/qPCR_data"
output_dir   <- "path/to/output_dir"
plate_width  <- 12
plate_height <- 8

plates_to_pooling(
    input_file,
    qpcr_dir,
    output_dir,
    plate_width,
    plate_height
)
```

### plates_to_pooling() arguments
- `input_file`: A path to the Excel file that was the output of the meta_to_plates() function.
- `qpcr_dir`: Directory with qPCR data. Filenames should be in a specific format mentioned below.
- `output_dir`: Directory where you would like all output files. The output files are explained in more detail below.
- `plate_width`: The number of columns per plate. Default = 12. This should be the same value you used in meta_to_plates().
- `plate_height`: The number of rows per plate. Default = 8. Max allowed = 13. This should be the same value you used in meta_to_plates().

### plates_to_pooling() input
- `input_file`: This is the output file from the meta_to_plates() function. An example can be viewed at `test_data/results.xlsx`
- `qpcr_dir`: An example of the qpcr directory can be viewed at `test_data/qPCR_test_data`
  - Each file in the qPCR directory should have names that end in `_$assay_$plate.txt`. An example name would be `20230901_Extraction_16S_Plate1.txt`. 
  - These txt files should be tab seperated files with the columns `Experiment_name`, `Position`, `Sample`, and `EPF`.

### plates_to_pooling() output
The output directory will contain three .csv files and multiple .pdf files.

- `reps_to_discard.csv`: This file has information on replicates that have been flaged for discarding.
  - `assay`: This column will have the assay of the discard replicate.
  - `plate_number`: The plate number of the discard replicate.
  - `pos`: The position on the plate of the discard replicate.
  - `sample_replicate`: The sample number and the replicate number of the discard replicate.
  - `discard`: The flag DISCARD.

- `biomek_pooling_workbook.csv` &
- `biomek_pooling_workbook_mamaaaaa.csv`: I currently don't understand the differernce between these files.

- `raw_$assay_$plate.pdf`: Plate plots before cleaning.

- `clean_$assay_$plate.pdf`: Plate plots after cleaning.
