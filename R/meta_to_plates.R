source("/data/sandbox/adam/qPCR_setup_pooling/R/import_index_df.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/import_meta_df.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/get_value_when_num_le_cutoff.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/get_possible_counts.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/reformat_meta_df.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/create_position_df.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/get_plate_dfs.R")
source("/data/sandbox/adam/qPCR_setup_pooling/R/export_plates_to_excel.R")
library(tidyverse)
library(readxl)
library(openxlsx)

meta_to_plates <- function(input_file,
                           output_file,
                           assays,
                           run,
                           plate_width = 12,
                           plate_height = 8,
                           controls = c("ITC", "NTC"),
                           control_pattern = "WC|DI|EB|BC|NTC|ITC|Cont") {

    meta_df  <- import_meta_df(input_file)
    index_df <- import_index_df(input_file, assays)

    # Get sample ids
    sample_ids <- meta_df$SAMPLEID

    # Make sure there is no duplicate samples
    if (length(sample_ids) != length(unique(sample_ids))) {
        duplicates <- meta_df[
            duplicated(meta_df$SAMPLE_ID) |
            duplicated(meta_df$SAMPLE_ID, fromLast = TRUE),
            "SAMPLE_ID"
        ]
        stop(
            paste0(
                "You can't have duplicate sample ids. Duplicates: ",
                duplicates
            )
        )
    }

    # Make sure the plate_height is a valid number
    if (plate_height > 13) {
        stop(
            paste0(
                "Plate height can't be greater than 13. Plate height: ",
                plate_height
            )
        )
    }

    # The number of plates and primers we need
    # is determined by the number of samples we have
    sample_count <- length(sample_ids)
    possible_counts <- get_possible_counts(
        sample_ids,
        controls,
        plate_height,
        plate_width
    )
    possible_plate_counts <- possible_counts$plates
    possible_fw_counts <- possible_counts$fws
    possible_rv_counts <- possible_counts$rvs
    sample_cutoffs <- possible_counts$cutoffs

    tryCatch(
        plate_count  <- get_value_when_num_le_cutoff(
            sample_count,
            sample_cutoffs,
            possible_plate_counts
        ),
        error = function(e){
            message("Error: You may have too many samples\n", e)
        }
    )

    fw_count <- get_value_when_num_le_cutoff(
            sample_count,
            sample_cutoffs,
            possible_fw_counts
        )

    rv_count <- get_value_when_num_le_cutoff(
            sample_count,
            sample_cutoffs,
            possible_rv_counts
        )

    meta_df     <- reformat_meta_df(
        meta_df,
        index_df,
        assays,
        controls,
        plate_count,
        fw_count,
        rv_count,
        plate_height,
        plate_width
    )

    print(meta_df$RV_NO)

    position_df <- create_position_df(
        meta_df,
        assays,
        plate_count,
        plate_height,
        plate_width,
        control_pattern
    )

    # plates and big_plates will be lists of data frames
    plates       <- list()
    big_plates   <- list()

    for (assay in assays) {
        first_fw     <- 1
        last_fw      <- plate_width
        first_rv     <- 1
        last_rv      <- plate_height
        first_sample <- 1
        last_sample  <- plate_height*plate_width
        for (plate_num in (1:plate_count)) {
            plate_dfs <- get_plate_dfs(
                first_fw,
                last_fw,
                first_rv,
                last_rv,
                first_sample,
                last_sample,
                sample_ids,
                index_df,
                fw_count,
                rv_count,
                assay,
                plate_num,
                plate_height,
                plate_width,
                meta_df
            )

            key                  <- paste0(assay, plate_num)
            plates[key][[1]]     <- plate_dfs$small_plate
            big_plates[key][[1]] <- plate_dfs$big_plate

            first_sample <- first_sample + (plate_height*plate_width)
            last_sample  <- last_sample + (plate_height*plate_width)

            if ((last_rv + plate_height) > rv_count) {
                first_rv    <- 1
                last_rv     <- plate_height
                first_fw    <- first_fw + plate_width
                last_fw     <- last_fw + plate_width

            } else {
                first_rv    <- first_rv + plate_height
                last_rv     <- last_rv + plate_height
            }
        }
    }

    export_plates_to_excel(
        assays,
        plates,
        big_plates,
        meta_df,
        position_df,
        output_file
    )
}