source("R/import_index_df.R")
source("R/import_meta_df.R")
source("R/get_value_when_num_le_cutoff.R")
source("R/get_possible_counts.R")
source("R/reformat_meta_df.R")
source("R/create_position_df.R")
source("R/get_plate_dfs.R")
source("R/export_plates_to_excel.R")
suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(openxlsx))

meta_to_plates <- function(input_file,
                           output_file,
                           assays,
                           run,
                           plate_width = 12,
                           plate_height = 8,
                           controls = c("ITC", "NTC"),
                           control_pattern = "WC|DI|EB|BC|NTC|ITC|Cont") {

    # Import our data
    meta_df  <- import_meta_df(input_file, run)
    index_df <- import_index_df(input_file, assays)

    # Remove invisible columns in Excel
    meta_df  <- meta_df %>% select(-contains("..."))
    # Remove control samples that shouldn't be in the input metadata
    for (control in controls) {
        meta_df <- meta_df[!grepl(control, meta_df$SAMPLEID),]
    }

    # Get sample ids
    sample_ids <- meta_df$SAMPLEID

    # Make sure the plate_height is a valid number
    if (plate_height > 13) {
        stop(
            paste0(
                "Plate height can't be greater than 13. Plate height: ",
                plate_height
            )
        )
    }

    # Get vectors to be used to get `plate_count`, `fw_count`, and `rv_count`
    # These vectors are:
    # - a vector of sample cutoffs based on based on number of samples per plate
    #   E.g., c(96, 192, 288, 384)
    # - a vector of possible plate counts
    #   E.g., c(1, 2, 3, 4)
    # - a vector of possible fw counts
    #   E.g., c(12, 12, 12, 24)
    # - a vector of possible rv counts
    #   E.g., c(8, 16, 24, 24)
    sample_count <- length(sample_ids)
    possible_counts <- get_possible_counts(
        sample_ids,
        controls,
        plate_height,
        plate_width
    )
    possible_plate_counts <- possible_counts$plates
    possible_fw_counts    <- possible_counts$fws
    possible_rv_counts    <- possible_counts$rvs
    sample_cutoffs        <- possible_counts$cutoffs

    # We should only need this tryCatch once
    tryCatch(
        plate_count  <- get_value_when_num_le_cutoff(
            sample_count,
            sample_cutoffs,
            possible_plate_counts
        ),
        error = function(e) {
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

    meta_df <- reformat_meta_df(
        meta_df,
        index_df,
        assays,
        controls,
        plate_count,
        fw_count,
        rv_count,
        plate_height,
        plate_width,
        run
    )

    position_df <- create_position_df(
        meta_df,
        assays,
        plate_count,
        plate_height,
        plate_width,
        control_pattern
    )

    # plates and big_plates will be named lists of data frames
    plates       <- list()
    big_plates   <- list()

    for (assay in assays) {
        first_fw     <- 1
        last_fw      <- plate_width
        first_rv     <- 1
        last_rv      <- plate_height
        first_sample <- 1
        last_sample  <- plate_height * plate_width
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

            first_sample <- first_sample + (plate_height * plate_width)
            last_sample  <- last_sample + (plate_height * plate_width)

            # We ran out of rv primers, so move onto next set of fw primers
            #
            # TODO: I may need to test this more to make sure it works
            #
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
        output_file,
        plate_height
    )
}
