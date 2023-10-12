######################################
# Setup
######################################

library("readxl")
library("openxlsx")
library("tidyverse")
library("dplyr")
library("gtools")


#### Define your input file name ####
# This should be in excel file format and contain a sheet titled" metadata"
# with columns called "sample_id" and "sequencing_run"
# There should also be a sheet for each assay you intend you use (e.g., 16S, MiFish etc)
# These index sheets should contain the name of the assay and the word index (e.g., 16S_index) 
# Each index sheet should have columns 'primer_#', 'primer_seq', 'tags', and 'fw_rv'
# the 'fw_rv' column should have 'fw' or 'rv' string values 
excel_file  <- "test_data_complexity.xlsx"
#### Define your output file name ####
output_file <- "test_data_output.xlsx"
#### Define which assays you intend to use #### 
assays      <- c("16S", 'MiFish')
#### Define which selection of samples in your metadata sheet you're selecting for this set-up ####
run         <- "run2"

######################################
# Function
######################################
# This function sets up the layout of the plates both 96 and 384
fill_plates_and_meta <- function(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df) {
    # Get a vector of just the primers we need
    fw_primers <- c(curr_index_df[curr_index_df$FWRV == "FW", "PRIMER#"][1:fw_count])
    rv_primers <- c(curr_index_df[curr_index_df$FWRV == "RV", "PRIMER#"][1:rv_count])
    
    if (length(fw_primers) < fw_count) {
        stop(paste0("You need at least ", fw_count, " fw primers. Only found ", length(fw_primers)))
    }
    if (length(rv_primers) < rv_count) {
        stop(paste0("You need at least ", rv_count, " rv primers. Only found ", length(rv_primers)))
    }

    plate     <- data.frame(matrix(ncol = 12, nrow = 8))
    big_plate <- data.frame(matrix(ncol = 24, nrow = 16))

    # We need to duplicate the primers for the big_plates, but we can't have duplicates yet, so let's append _1, or _2
    doubled_fw_primers <- rep(fw_primers[first_fw:last_fw], each = 2)
    doubled_fw_primers <- ifelse(seq_along(doubled_fw_primers) %% 2 == 1, paste0(doubled_fw_primers, "_1"), paste0(doubled_fw_primers, "_2"))
    doubled_rv_primers <- rep(rv_primers[first_rv:last_rv], each = 2)
    doubled_rv_primers <- ifelse(seq_along(doubled_rv_primers) %% 2 == 1, paste0(doubled_rv_primers, "_1"), paste0(doubled_rv_primers, "_2"))
            
    colnames(plate)     <- fw_primers[first_fw:last_fw]
    colnames(big_plate) <- doubled_fw_primers
    rownames(plate)     <- rv_primers[first_rv:last_rv]
    rownames(big_plate) <- doubled_rv_primers

    # Track the well we are up to so we can include the well info
    well_rows      <- c("A", "B", "C", "D", "E", "F", "G", "H")
    well_col       <- 0

    # There is 12 columns and 8 rows per plate, 
    # so these for loops will iterate 96 times
    curr_sam_index <- first_sample
    for (col in colnames(plate)) {
        well_row_index <- 0
        well_col       <- well_col + 1
        for (row in rownames(plate)) {
            well_row_index <- well_row_index + 1
            well           <- paste0(well_rows[well_row_index], well_col)
            if (curr_sam_index <= last_sample) {
                curr_sample     <- sample_ids[curr_sam_index]

                if (!is.na(curr_sample)) {
                    # Add the current sample to the current cell
                    plate[row,col] <- curr_sample

                    # Each sample should be added to the big plates three times
                    big_plate[paste0(row, "_1"),paste0(col, "_1")] <- paste0(curr_sample, "-1")
                    big_plate[paste0(row, "_1"),paste0(col, "_2")] <- paste0(curr_sample, "-2")
                    big_plate[paste0(row, "_2"),paste0(col, "_1")] <- paste0(curr_sample, "-3")
                    big_plate[paste0(row, "_2"),paste0(col, "_2")] <- paste0(curr_sample, "-pool")

                    # Add primer info to metadata
                    meta_row <- which(meta_df$SAMPLEID == curr_sample & meta_df$ASSAY == assay)
                    meta_df[meta_row, "FW_NO"]       <- col
                    meta_df[meta_row, "RV_NO"]       <- row
                    meta_df[meta_row, "FW_TAG"]      <- subset(curr_index_df, `PRIMER#` == col & FWRV == "FW", select = TAGS)$TAGS
                    meta_df[meta_row, "RV_TAG"]      <- subset(curr_index_df, `PRIMER#` == row & FWRV == "RV", select = TAGS)$TAGS
                    meta_df[meta_row, "FW_FULL_SEQ"] <- subset(curr_index_df, `PRIMER#` == col & FWRV == "FW", select = PRIMERSEQ)$PRIMERSEQ
                    meta_df[meta_row, "RV_FULL_SEQ"] <- subset(curr_index_df, `PRIMER#` == row & FWRV == "RV", select = PRIMERSEQ)$PRIMERSEQ
                    meta_df[meta_row, "PLATE"]       <- plate_num
                    meta_df[meta_row, "WELL"]        <- well
                }
                        
                curr_sam_index <- curr_sam_index + 1
            }
        }
    }

    # Track the well and replicate we are up to \
    well_rows      <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
    well_row_index <- 1
    well_cols      <- 1:24
    well_col_index <- 1
    replicate      <- "1"

    # There is 24 columns and 16 rows per big plate, 
    # so we need to loop 384 times
    curr_sam_index <- first_sample
    for (i in 1:384) {
        well <- paste0(well_rows[well_row_index], well_cols[well_col_index])
        if (curr_sam_index <= last_sample) {
            curr_sample     <- sample_ids[curr_sam_index]

            # Samples with these patterns should be flaged as control samples
            if (grepl("WC|DI|EB|BC|NTC|ITC|Cont", c(curr_sample))) {
                sample_type <- "control"
            } else {
                sample_type <- "sample"
            }

            # Add position info to position df
            pos_row <- which(position_df$sample_id == curr_sample & position_df$replicate == replicate & position_df$ASSAY == assay)
            position_df[pos_row, "Pos"]              <- well
            position_df[pos_row, "sample_replicate"] <- paste0(curr_sample, "-", replicate)
            position_df[pos_row, "plate_number"]     <- paste0("Plate", plate_num)
            position_df[pos_row, "assay_plate_number_Pos"] <- paste0(assay, ".", "Plate", plate_num, ".", well)
            position_df[pos_row, "plate_number_Pos"] <- paste0("Plate", plate_num, ".", well)
            position_df[pos_row, "sample_type"]      <- sample_type

            # If the current replicate is '1', we need to move one column to the right and stay in the same row
            if (replicate == "1") {
                replicate      <- "2"
                well_col_index <- well_col_index + 1

            # If the current replicate is '2', we need to move one column to the left and one row down
            } else if (replicate == "2") {
                replicate      <- "3"
                well_row_index <- well_row_index + 1
                well_col_index <- well_col_index - 1

            # If the current replicate is '3', we need to move one column to the right
            } else if (replicate == "3") {
                replicate      <- "pool"
                well_col_index <- well_col_index + 1  
                
                # If the current replicate is 'pool', we need to check if we are up to the last row 
            } else if (replicate == "pool") {
              replicate      <- "1"
              curr_sam_index <- curr_sam_index + 1 

                 # If it is the last row, move to the first row, and move across one columns
                if (well_row_index == 16) {
                    well_row_index <- 1
                    well_col_index <- well_col_index + 1

                # If it is not the last row, move down one row, but stay in the same column    
                } else {
                    well_row_index <- well_row_index + 1
                    well_col_index <- well_col_index - 1
                    
                }
            }
        }
    }

    list(df1 = plate, df2 = big_plate, df3 = meta_df, df4 = position_df)
}

######################################
# Main
######################################
# Validates and inputs the data from your input file to the layout defined above
# Get sheet named 'metadata'
sheets          <- excel_sheets(excel_file)
curr_index_name <- grep(glob2rx('metadata'), sheets, ignore.case = TRUE, value = TRUE)

# Validate number of sheets called 'metadata'
if (length(curr_index_name) == 0) {
    stop("No sheets found nemed 'metadata'")
} 
if (length(curr_index_name) > 1) {
    stop("Multiple sheets found named 'metadata'")
}

# Import metadata
meta_df   <- read_excel(excel_file, sheet="metadata") %>%
    as.data.frame()
names(meta_df) <- gsub(" ", "", names(meta_df)) %>%
    { gsub("_", "", .) } %>%
    toupper()

# Add empty columns
meta_df["FW_NO"]       <- NULL
meta_df["RV_NO"]       <- NULL
meta_df["FW_TAG"]      <- NULL
meta_df["RV_TAG"]      <- NULL
meta_df["FW_FULL_SEQ"] <- NULL
meta_df["RV_FULL_SEQ"] <- NULL
meta_df["ASSAY"]       <- NULL
meta_df["PLATE"]       <- NULL
meta_df["WELL"]        <- NULL

# Validata metadata columns
if (!("SAMPLEID" %in% colnames(meta_df))) {
    stop(paste0("Sheet 'metadata' must contain a 'sample_id' column"))
}

# Subset data if there is a sequencing run column
if ("SEQUENCINGRUN" %in% colnames(meta_df)) {
    meta_df <- subset(meta_df, SEQUENCINGRUN == run)
}

# Get sample ids
sample_ids             <- meta_df$SAMPLEID

# Make sure there is no duplicate samples
if (length(sample_ids) != length(unique(sample_ids))) {
    duplicates <- meta_df[duplicated(meta_df$SAMPLE_ID) | duplicated(meta_df$SAMPLE_ID, fromLast = TRUE), "SAMPLE_ID"]
    stop(paste0("You can't have duplicate sample ids. Duplicates: ", duplicates))
}

# plates and big_plates will be lists of data frames
plates       <- list()
big_plates   <- list()

# We need to duplicate the samples in the metadata for each assay
meta_df            <- meta_df[rep(seq_len(nrow(meta_df)), each = length(assays)), ]
meta_df$ASSAY      <- rep(assays, length.out = nrow(meta_df))
row.names(meta_df) <- NULL

# Create a position_df to hold information that can link the samples to the big plate positions
position_df                     <- meta_df[, c("SAMPLEID", "ASSAY")]
colnames(position_df)[colnames(position_df) == 'SAMPLEID'] <- 'sample_id'
position_df["Pos"]              <- NULL
position_df["sample_replicate"] <- NULL
position_df["replicate"]        <- NULL
position_df["plate_number"]     <- NULL
position_df["assay_plate_number_Pos"] <- NULL
position_df["plate_number_Pos"] <- NULL
position_df["sample_type"]      <- NULL
replicates                      <- c("1", "2", "3", "pool")
position_df                     <- position_df[rep(seq_len(nrow(position_df)), each = length(replicates)), ]
position_df$replicate           <- rep(replicates, length.out = nrow(position_df))
row.names(position_df)          <- NULL

for (assay in assays) {
    # Get sheet names that match the correct pattern
    curr_index_sheet <- grep(glob2rx(paste0('*', assay, '*index*|*index*', assay, '*')), sheets, ignore.case = TRUE, value = TRUE)

    # Validate number of sheets that matched pattern
    if (length(curr_index_sheet) == 0) {
        stop(paste0("No sheets found containing words 'index' and ", assay))
    } 
    if (length(curr_index_sheet) > 1) {
        stop(paste0("Multiple sheets found containing words 'index' and ", assay))
    }

    # Import index sheet
    curr_index_df <- read_excel(excel_file, sheet=curr_index_sheet) %>%
        as.data.frame()
    names(curr_index_df) <- gsub(" ", "", names(curr_index_df)) %>%
        { gsub("_", "", .) } %>%
        toupper()
    
    # Validate index sheet columns
    if (!("PRIMER#" %in% colnames(curr_index_df))) {
      stop(paste0("Sheet '", curr_index_sheet, "' must contain a 'primer_#' column"))
    }
    # Validate index sheet columns
    if (!("PRIMERSEQ" %in% colnames(curr_index_df))) {
      stop(paste0("Sheet '", curr_index_sheet, "' must contain a 'primer_seq' column"))
    }
    # Validate index sheet columns
    if (!("TAGS" %in% colnames(curr_index_df))) {
      stop(paste0("Sheet '", curr_index_sheet, "' must contain a 'tags' column"))
    }
    # Validate index sheet columns
    if (!("FWRV" %in% colnames(curr_index_df))) {
      stop(paste0("Sheet '", curr_index_sheet, "' must contain a 'fw_rv' column"))
    }

    # Make sure there is no duplicate primers
    if (length(curr_index_df$`PRIMER#`) != length(unique(curr_index_df$`PRIMER#`))) {
        duplicates <- curr_index_df[duplicated(curr_index_df$`PRIMER#`) | duplicated(curr_index_df$`PRIMER#`, fromLast = TRUE), "PRIMER_#"]
        stop(paste0("You can't have duplicate primers. Duplicates: ", duplicates))
    }
    if (length(curr_index_df$PRIMERSEQ) != length(unique(curr_index_df$PRIMERSEQ))) {
        duplicates <- curr_index_df[duplicated(curr_index_df$PRIMERSEQ) | duplicated(curr_index_df$PRIMERSEQ, fromLast = TRUE), "PRIMER_SEQ"]
        stop(paste0("You can't have duplicate primer sequences. Duplicates: ", duplicates))
    }
    if (length(curr_index_df$TAGS) != length(unique(curr_index_df$TAGS))) {
        duplicates <- curr_index_df[duplicated(curr_index_df$TAGS) | duplicated(curr_index_df$TAGS, fromLast = TRUE), "TAGS"]
        stop(paste0("You can't have duplicate tags. Duplicates: ", duplicates))
    }

    # Make sure fw_rv column only has 'fw' or 'rv' values
    curr_index_df$FWRV <- toupper(curr_index_df$FWRV)
    unique_fw_rv       <- unique(curr_index_df$FWRV)
    for (i in unique_fw_rv) {
        if (i != "FW" & i != "RV") {
            stop(paste0("You can't have ", i, " in fw_rv column. Valid values are 'fw' or 'rv'"))
        }
    }

    # The number of plates and primers we need is determined by the number of samples we have
    if (length(sample_ids) <= 94) {
        plate_count <- 1
        fw_count    <- 12
        rv_count    <- 8
    } else if (length(sample_ids) <= 188) {
        plate_count <- 2
        fw_count    <- 12
        rv_count    <- 16
    } else if (length(sample_ids) <= 282) {
        plate_count <- 3
        fw_count    <- 12
        rv_count    <- 24
    } else if (length(sample_ids) <= 376) {
        plate_count <- 4
        fw_count    <- 24
        rv_count    <- 24
    } else if (length(sample_ids) <= 470) {
        plate_count <- 5
        fw_count    <- 24
        rv_count    <- 24
    } else if (length(sample_ids) <= 564) {
        plate_count <- 6
        fw_count    <- 24
        rv_count    <- 24
    } else if (length(sample_ids) <= 658) {
        plate_count <- 7
        fw_count    <- 36
        rv_count    <- 24
    } else if (length(sample_ids) <= 752) {
        plate_count <- 8
        fw_count    <- 36
        rv_count    <- 24
    } else if (length(sample_ids) <= 846) {
        plate_count <- 9
        fw_count    <- 36
        rv_count    <- 24
    } else {
        stop(paste0("You have ", length(sample_ids), " samples. Max allowed is 8"))
    }

    # Loop through the plate numbers and use that number with the assay as a key to store our data
    for (plate_num in (1:plate_count)) {
        key <- paste0(assay, plate_num)
        # For the first plate of the assay, we start at the first sample and end at sample 95 because we want to
        # leave the bottom right corner empty
        if (plate_num == 1) {
            first_fw     <- 1
            last_fw      <- 12
            first_rv     <- 1
            last_rv      <- 8
            first_sample <- 1
            last_sample  <- 94
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
            
        # for the other plates, we are just changing the numbers a bit
        } else if (plate_num == 2) {
            first_fw     <- 1
            last_fw      <- 12
            first_rv     <- 9
            last_rv      <- 16
            first_sample <- 95
            last_sample  <- 188
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
            
        } else if (plate_num == 3) {
            first_fw     <- 1
            last_fw      <- 12
            first_rv     <- 17
            last_rv      <- 24
            first_sample <- 189
            last_sample  <- 283
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4

        } else if (plate_num == 4) {
            first_fw     <- 13
            last_fw      <- 24
            first_rv     <- 1
            last_rv      <- 8
            first_sample <- 284
            last_sample  <- 378
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        } else if (plate_num == 5) {
            first_fw     <- 13
            last_fw      <- 24
            first_rv     <- 9
            last_rv      <- 16
            first_sample <- 379
            last_sample  <- 473
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        } else if (plate_num == 6) {
            first_fw     <- 25
            last_fw      <- 36
            first_rv     <- 17
            last_rv      <- 24
            first_sample <- 474
            last_sample  <- 567
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        } else if (plate_num == 7) {
            first_fw     <- 25
            last_fw      <- 36
            first_rv     <- 1
            last_rv      <- 8
            first_sample <- 568
            last_sample  <- 661
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        } else if (plate_num == 8) {
            first_fw     <- 25
            last_fw      <- 36
            first_rv     <- 9
            last_rv      <- 16
            first_sample <- 662
            last_sample  <- 757
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        } else if (plate_num == 9) {
            first_fw     <- 13
            last_fw      <- 24
            first_rv     <- 17
            last_rv      <- 24
            first_sample <- 758
            last_sample  <- 851
            dfs          <- fill_plates_and_meta(first_fw, last_fw, first_rv, last_rv, first_sample, last_sample, meta_df, sample_ids, curr_index_df, fw_count, rv_count, assay, plate_num, position_df)
            plates[key][[1]]     <- dfs$df1
            big_plates[key][[1]] <- dfs$df2
            meta_df              <- dfs$df3
            position_df          <- dfs$df4
        }
    }
}

# Create a workbook object
wb <- createWorkbook()

# Add a sheet to the workbook
addWorksheet(wb, "plates")
addWorksheet(wb, "big_plates")
addWorksheet(wb, "metadata")
addWorksheet(wb, "position_df")

# We will loop through the assays and plate numbers and use those values as keys
# to get the correct data for populating the output Excel file
# 'row_start' and 'big_row_start' are to track the rows we are up to in the Excel file
# so that the next plate is placed under the previous plate in the Excel sheet
row_start        <- 1
big_row_start    <- 1

for (assay in assays) {
    for (plate_num in 1:plate_count) {
        key <- paste0(assay, plate_num)

        # Move the row names and col names to the actual df, remove NAs, 
        # and remove the _1 and _2 we appended to some of those names earlier
        plates[key][[1]]                   <- rownames_to_column(plates[key][[1]], var = paste0(assay, " plate ", plate_num))
        big_plates[key][[1]]               <- rownames_to_column(big_plates[key][[1]], var = paste0(assay, " plate ", plate_num))
        big_plates[key][[1]][, 1]          <- gsub("_[12]$", "", big_plates[key][[1]][, 1])
        colnames(big_plates[key][[1]])[-1] <- gsub("_[12]$", "", colnames(big_plates[key][[1]])[-1])
        plates[key][[1]]                   <- rbind(colnames(plates[key][[1]]), plates[key][[1]])
        colnames(plates[key][[1]])         <- ""
        big_plates[key][[1]]               <- rbind(colnames(big_plates[key][[1]]), big_plates[key][[1]])
        colnames(big_plates[key][[1]])     <- ""
        plates[key][[1]]                   <- replace(plates[key][[1]], is.na(plates[key][[1]]) | plates[key][[1]] == "NA", "")
        big_plates[key][[1]]               <- replace(big_plates[key][[1]], is.na(big_plates[key][[1]]) | big_plates[key][[1]] == "NA", "")

        # Write the data to the sheet
        writeData(wb, sheet = "plates", x = plates[key][[1]], startRow = row_start, startCol = 1, colNames = FALSE, rowNames = FALSE)
        writeData(wb, sheet = "big_plates", x = big_plates[key][[1]], startRow = big_row_start, startCol = 1, colNames = FALSE, rowNames = FALSE)
        row_start     <- row_start + 10
        big_row_start <- big_row_start + 18
    }
}

# We can't forget the add the metadata which now has the primer info
writeData(wb, sheet = "metadata", meta_df, startRow = 1, startCol = 1)
writeData(wb, sheet = "position_df", position_df, startRow = 1, startCol = 1)
    
# Save the workbook as an Excel file
saveWorkbook(wb, output_file, overwrite = TRUE)

print(paste0("Finished! Output file: ", output_file))

# Calculate number of reactions per unique Fw and Rv primers
#### Forward primers ####
n_rx_unique_fw <- meta_df %>%
  group_by(ASSAY, FW_NO) %>%
  dplyr::summarise(n = n()) %>% # number of samples
  mutate(number_rxn = n * 3) %>% # number of samples in triplicate - rxn for NTC not included in this calculation
  print(n = Inf)

summary(n_rx_unique_fw$number_rxn)
# 24 or 72 rxn per primer (8 or 24 samples in triplicate, respectively)

n_rx_unique_fw_single_plate <- n_rx_unique_fw %>%
  dplyr::filter(n <=8) %>% # single plate = 8 samples max in total
  print(n = Inf)

n_rx_unique_fw_single_plate_16s <- n_rx_unique_fw_single_plate %>%
  dplyr::filter(ASSAY == "16S")

n_rx_unique_fw_single_plate_16s$FW_NO[mixedorder(n_rx_unique_fw_single_plate_16s$FW_NO)]

n_rx_unique_fw_single_plate_COI <- n_rx_unique_fw_single_plate %>%
  dplyr::filter(ASSAY == "COI")

n_rx_unique_fw_single_plate_COI$FW_NO[mixedorder(n_rx_unique_fw_single_plate_COI$FW_NO)]

n_rx_unique_fw_multiple_plate <- n_rx_unique_fw %>%
  dplyr::filter(n > 8) %>% # multiple plates > 9 samples in total
  print(n = Inf)

n_rx_unique_fw_multiple_plate_16s <- n_rx_unique_fw_multiple_plate %>%
  dplyr::filter(ASSAY == "16S")

n_rx_unique_fw_multiple_plate_16s$FW_NO[mixedorder(n_rx_unique_fw_multiple_plate_16s$FW_NO)]

n_rx_unique_fw_multiple_plate_COI <- n_rx_unique_fw_multiple_plate %>%
  dplyr::filter(ASSAY == "COI")

n_rx_unique_fw_multiple_plate_COI$FW_NO[mixedorder(n_rx_unique_fw_multiple_plate_COI$FW_NO)]

#### Reverse primers ####
n_rx_unique_rv <- meta_df %>%
  group_by(ASSAY, RV_NO) %>%
  dplyr::summarise(n = n()) %>% # number of samples
  mutate(number_rxn = n * 3) %>% # number of samples in triplicate - rxn for NTC not included in this calculation
  print(n = Inf)

summary(n_rx_unique_rv$number_rxn)
# 36 or 72 rxn per primer (12 or 24 samples in triplicate)

n_rx_unique_rv_single_plate <- n_rx_unique_rv %>%
  dplyr::filter(n <=12) %>% # single plate = 12 samples max in total
  print(n = Inf)

n_rx_unique_rv_single_plate_16s <- n_rx_unique_rv_single_plate %>%
  dplyr::filter(ASSAY == "16S")

n_rx_unique_rv_single_plate_16s$RV_NO[mixedorder(n_rx_unique_rv_single_plate_16s$RV_NO)]

n_rx_unique_rv_single_plate_COI <- n_rx_unique_rv_single_plate %>%
  dplyr::filter(ASSAY == "COI")

n_rx_unique_rv_single_plate_COI$RV_NO[mixedorder(n_rx_unique_rv_single_plate_COI$RV_NO)]

n_rx_unique_rv_multiple_plate <- n_rx_unique_rv %>%
  dplyr::filter(n > 12) %>% # multiple plates > 13 samples in total
  print(n = Inf)

n_rx_unique_rv_multiple_plate_16s <- n_rx_unique_rv_multiple_plate %>%
  dplyr::filter(ASSAY == "16S")

n_rx_unique_rv_multiple_plate_16s$RV_NO[mixedorder(n_rx_unique_rv_multiple_plate_16s$RV_NO)]

n_rx_unique_rv_multiple_plate_COI <- n_rx_unique_rv_multiple_plate %>%
  dplyr::filter(ASSAY == "COI")

n_rx_unique_rv_multiple_plate_COI$RV_NO[mixedorder(n_rx_unique_rv_multiple_plate_COI$RV_NO)]

