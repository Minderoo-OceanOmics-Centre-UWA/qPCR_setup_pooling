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

meta_to_plates <- function(metadata,
                           index_file,
                           output_file,
                           assays,
                           run = "run1",
                           plate_width = 12,
                           plate_height = 8,
                           controls = c("NTC", "ITC"),
                           control_pattern = "WC|DI|EB|BC|NTC|ITC|Cont|BL",
                           skip_plates = list(),
                           skip_samples = list(),
                           strategy = "UDI") {
  
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
    if (length(metadata) == 1) {
        meta_df <- import_meta_df(metadata, run)
    } else {
        meta_df <- data.frame() 
        for (file in metadata) {
            tmp_df  <- import_meta_df(file, run)
            meta_df <- rbind(meta_df, tmp_df)
        }
    }
    index_df <- import_index_df(index_file, assays)
  
    # Remove invisible columns in Excel
    meta_df  <- meta_df %>% dplyr::select(-contains("..."))
    # Remove control samples that shouldn't be in the input metadata
    for (control in controls) {
        meta_df <- data.frame(meta_df[!grepl(control, meta_df$samp_name),])
        if (length(colnames(meta_df)) == 1) {
            colnames(meta_df) <- "samp_name"
        }
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
            fake_meta_df <- data.frame(samp_name = sample_ids[assay][[1]])
            tmp_meta_df  <- rbind.fill(fake_meta_df, tmp_meta_df)
        }
      
        sample_ids[assay][[1]] <- tmp_meta_df$samp_name
      
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
        sample_ids,
        strategy
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
      sample      <- meta_df[row, "samp_name"]
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
      filter(!grepl("^FAKESAMPLE_", samp_name))
    
    position_df <- position_df %>%
      filter(!grepl("^FAKESAMPLE_", sample_id))


  num_rows <- 384
    well_chars <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
    well_nums <- c(1:24)
    well_positions <- c()
    for (char in well_chars) {
        for (num in well_nums) {
            well_positions <- c(well_positions, paste0(char, num))
        }
    }
    
    for (assay in assays) {
        QS7_outfile <- paste0("QS7_plate_import_", assay, ".csv")
        
        write("* Block Type = 384-Well Block", file = QS7_outfile)
        write(paste0("* Date Created = ", format(Sys.time(), "%a %b %d %H:%M:%S AWST %Y")), file = QS7_outfile, append = TRUE)
        write("* Passive Reference = ", file = QS7_outfile, append = TRUE)
        write("* Barcode = 1237", file = QS7_outfile, append = TRUE)
        write("", file = QS7_outfile, append = TRUE)
        write("[Sample Setup]", file = QS7_outfile, append = TRUE)
        
        position_df_QS7 <- data.frame(
            "Well" = numeric(num_rows),
            "Well Position" = character(num_rows),
            "Sample Name" = character(num_rows),
            "Sample Color" = character(num_rows),
            "Biogroup Name" = character(num_rows),
            "Biogroup Color" = character(num_rows),
            "Target Name" = character(num_rows),
            "Target Color" = character(num_rows),
            "Task" = character(num_rows),
            "Reporter" = character(num_rows),
            "Quencher" = character(num_rows),
            "Quality" = character(num_rows),
            "Quantity" = character(num_rows),
            "Comments" = character(num_rows),
            check.names = FALSE 
        )
        
        position_df_QS7["Well"] <- c(1:384)
        position_df_QS7["Well Position"] <- well_positions
        
        filled_positions <- position_df$Pos
        
        for (pos in filled_positions) {
            row_QS7 = trimws(position_df_QS7[position_df_QS7["Well Position"] == pos][1])
            row_pdf = rownames(position_df[position_df$Pos == pos & position_df$assay == assay, ])
        
            position_df_QS7[row_QS7, "Sample Name"] <- position_df[row_pdf, "sample_id"]
            position_df_QS7[row_QS7, "Target Name"] <- assay
            position_df_QS7[row_QS7, "Target Color"] <- "RGB(86,214,243)"
            position_df_QS7[row_QS7, "Task"] <- "UNKNOWN"
            position_df_QS7[row_QS7, "Reporter"] <- "SYBR"
        }
        write.table(position_df_QS7, QS7_outfile, 
                    append = TRUE, 
                    sep = ",", 
                    row.names = FALSE, 
                    col.names = TRUE)
    }
  
    export_plates_to_excel(
        assays,
        plates,
        big_plates,
        meta_df,
        position_df,
        output_file,
        plate_height,
        plate_count,
        strategy
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
