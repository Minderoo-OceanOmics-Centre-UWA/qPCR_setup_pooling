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
suppressMessages(library(plyr))

meta_to_plates <- function(input_file,
                           output_file,
                           assays,
                           run,
                           plate_width = 12,
                           plate_height = 8,
                           controls = c("NTC", "ITC"),
                           control_pattern = "WC|DI|EB|BC|NTC|ITC|Cont|BL",
                           skip_plates = list(),
                           skip_samples = list()) {
  
    # Make sure the plate_height is a valid number
    if (plate_height > 13) {
        stop(
            paste0(
                "Plate height can't be greater than 13. Plate height: ",
                plate_height
            )
        )
    }
    
    # Import our data
    meta_df  <- import_meta_df(input_file, run)
    index_df <- import_index_df(input_file, assays)
  
    # Remove invisible columns in Excel
    meta_df  <- meta_df %>% dplyr::select(-contains("..."))
    # Remove control samples that shouldn't be in the input metadata
    for (control in controls) {
        meta_df <- meta_df[!grepl(control, meta_df$SAMPLEID),]
    }
    
    sample_ids            <- list()
    sample_count          <- list()
    possible_plate_counts <- list()
    possible_fw_counts    <- list()
    possible_rv_counts    <- list()
    sample_cutoffs        <- list()
    plate_count           <- list()
    fw_count              <- list()
    rv_count              <- list()
    skip_samples_param    <- skip_samples
  
    for (assay in assays) {
        tmp_meta_df <- meta_df
        sample_ids[assay] <- list(c())
      
        # Let's add fake samples for when the user wants to skip samples/plates
        if (skip_plates[assay][[1]] > 0 | skip_samples_param[assay][[1]] > 0) {
        
            skip_samples[assay][[1]] <- skip_samples_param[assay][[1]] + (((plate_height * plate_width) - length(controls)) * skip_plates[assay][[1]])
        
            for (i in 1:skip_samples[assay][[1]]) {
                sample_ids[assay][[1]] <- c(sample_ids[assay][[1]], paste0("FAKESAMPLE_", i)) 
            }
            fake_meta_df <- data.frame(SAMPLEID = sample_ids[assay][[1]], SEQUENCINGRUN = run)
            tmp_meta_df  <- rbind.fill(fake_meta_df, tmp_meta_df)
        }
      
        sample_ids[assay][[1]] <- tmp_meta_df$SAMPLEID
      
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
        sample_count[assay][[1]]    <- length(sample_ids[assay][[1]])
        possible_counts <- get_possible_counts(
            sample_ids[assay][[1]],
            controls,
            plate_height,
            plate_width
        )
        possible_plate_counts[assay][[1]] <- possible_counts$plates
        possible_fw_counts[assay][[1]]    <- possible_counts$fws
        possible_rv_counts[assay][[1]]    <- possible_counts$rvs
        sample_cutoffs[assay][[1]]        <- possible_counts$cutoffs
      
        # We should only need this tryCatch once
        tryCatch(
            plate_count[assay][[1]]  <- get_value_when_num_le_cutoff(
                sample_count[assay][[1]],
                sample_cutoffs[assay][[1]],
                possible_plate_counts[assay][[1]]
            ),
            error = function(e) {
                message("Error: You may have too many samples\n", e)
            }
        )
      
        fw_count[assay][[1]] <- get_value_when_num_le_cutoff(
            sample_count[assay][[1]],
            sample_cutoffs[assay][[1]],
            possible_fw_counts[assay][[1]]
        )
      
        rv_count[assay][[1]] <- get_value_when_num_le_cutoff(
            sample_count[assay][[1]],
            sample_cutoffs[assay][[1]],
            possible_rv_counts[assay][[1]]
        )
    }

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
        run,
        sample_ids
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
        for (plate_num in (1:plate_count[assay][[1]])) {
            plate_dfs <- get_plate_dfs(
                first_fw,
                last_fw,
                first_rv,
                last_rv,
                first_sample,
                last_sample,
                sample_ids[assay][[1]],
                index_df,
                fw_count[assay][[1]],
                rv_count[assay][[1]],
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
            if ((last_rv + plate_height) > rv_count[assay][[1]]) {
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
    
    meta_df$sample_type <- NA
    control_vect        <- str_split_1(control_pattern, "\\|")
    for (row in 1:nrow(meta_df)) {
      sample      <- meta_df[row, "SAMPLEID"]
      sample_type <- "sample"
      
      for (i in control_vect) {
        if (grepl(i, sample, fixed = TRUE)) {
          sample_type <- paste0(i, "_Control")
          break
        }
      }
      
      meta_df[row, "sample_type"] <- sample_type
    }
    
    meta_df <- meta_df %>%
      filter(!grepl("^FAKESAMPLE_", SAMPLEID))
    
    position_df <- position_df %>%
      filter(!grepl("^FAKESAMPLE_", sample_id))

    export_plates_to_excel(
        assays,
        plates,
        big_plates,
        meta_df,
        position_df,
        output_file,
        plate_height,
        plate_count
    )

    for (assay in assays) {
      last_plate <- plate_count[assay]
      biomek_out_csv <- data.frame(row.names=row.names(position_df[position_df$assay == assay & position_df$plate_number == paste0("Plate", last_plate) & position_df$replicate != "pool",]))
      biomek_out_csv$SourcePosition <- "Reservor"
      biomek_out_csv$Quarter <- 1
      biomek_out_csv$Volume <- 12.8
      biomek_out_csv$DestinationPosition <- "qPCRplate"
      biomek_out_csv$Dest_1 <- position_df[position_df$assay == assay & position_df$plate_number == paste0("Plate", last_plate) & position_df$replicate != "pool",]$Pos
      biomek_out_csv$ID <- ""
      
      write.csv(biomek_out_csv, paste0(assay, "_biomek_MM_Plating.csv"), row.names = FALSE)
    }
}
