reformat_meta_df <- function(meta_df,
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
                             strategy) {

    control_count <- length(controls)
    sample_count_per_plate <- (plate_height * plate_width) - control_count
    sample_rows <- list()
    sample_count <- list()
    
    for (assay in assays) {
        sample_rows[assay][[1]] <- data.frame(
            samp_name = sample_ids[assay][[1]]
        )
        sample_rows[assay][[1]]$TMP <- ""
        
        # Add control samples to meta_df.
        # We need to add the control samples at the end of each plate.
        left_of_controls <- sample_count_per_plate
        fake_plates <- 0
        for (plate_num in 1:plate_count[assay][[1]]) {
          if (plate_num == 1) {
              if (all(grepl("^FAKESAMPLE_", sample_rows[assay][[1]][1:left_of_controls, ]$samp_name))) {
                  fake_plates <- fake_plates + 1
              }  
          } else {
            if (all(grepl("^FAKESAMPLE_", sample_rows[assay][[1]][(prev_left_of_controls + control_count + 1):left_of_controls, ]$samp_name))) {
              fake_plates <- fake_plates + 1
            }  
          }
          
          if ((plate_num - fake_plates) == 0) {
            control_rows <- data.frame(
                samp_name = paste0("FAKESAMPLE_", controls, "_", plate_num)
            )
          } else {
            control_rows <- data.frame(
                samp_name = paste0(controls, "_", (plate_num - fake_plates))
            ) 
          }
          
          # This will be true for the last plate when there isn't enough
          # samples to fill the plate
          if (left_of_controls > nrow(sample_rows[assay][[1]])) {
            prev_left_of_controls <- left_of_controls
            left_of_controls <- nrow(sample_rows[assay][[1]])
            sample_rows[assay][[1]] <- rbind.fill(sample_rows[assay][[1]], control_rows)
            
            left_of_controls <- (
              left_of_controls +
                control_count +
                sample_count_per_plate
            )
            
            # This will be true for most plates
          } else {
            # We need to handle the last plate differently if the plate is filled completely with samples  
            if (plate_num == plate_count[assay][[1]] & left_of_controls == nrow(sample_rows[assay][[1]])) {
              sample_rows[assay][[1]] <- rbind.fill(
                sample_rows[assay][[1]][1:left_of_controls, ],
                control_rows
              )
                
              prev_left_of_controls <- left_of_controls
              left_of_controls <- (
                left_of_controls +
                control_count +
                sample_count_per_plate
              )
            } else {
              sample_rows[assay][[1]] <- rbind.fill(
                sample_rows[assay][[1]][1:left_of_controls, ],
                control_rows,
                sample_rows[assay][[1]][(left_of_controls + 1):nrow(sample_rows[assay][[1]]), ]
              )
                
              prev_left_of_controls <- left_of_controls
              left_of_controls <- (
                left_of_controls +
                control_count +
                sample_count_per_plate
              )
            }
          }
        }
        sample_ids[assay][[1]] <- sample_rows[assay][[1]]$samp_name
        sample_count[assay][[1]] <- length(sample_ids[assay][[1]])
        sample_rows[assay][[1]]$TMP <- NULL
    }

    

    # We need to duplicate the samples for each assay
    sample_rows <- rbind.fill(sample_rows)
    assay_vec   <- c()
    for (assay in assays) {
        assay_vec <- c(assay_vec, rep(assay, length.out = sample_count[assay][[1]]))
    }

    sample_rows$assay      <- assay_vec
    row.names(sample_rows) <- NULL
    
    # Create a vector of characters based on plate height.
    well_rows      <- c()
    char           <- "A"
    for (i in 1:plate_height) {
      well_rows <- append(well_rows, char)
      # Get next letter
      char      <- LETTERS[which(LETTERS == char) + 1L]
    }
    
    meta_df <- merge(sample_rows, meta_df, by = c("samp_name"), all.x = TRUE)

    meta_df <- meta_df[match(paste(sample_rows$samp_name, sample_rows$assay), 
                                    paste(meta_df$samp_name, meta_df$assay)), ]
    rownames(meta_df) <- NULL
    
    
    #tmp_meta_df <- meta_df
    
    
    #meta_df <- tmp_meta_df
    if (strategy == "UC") {
        for (curr_assay in assays) {
            # Get vectors of the primers
            fw_primers <- c(index_df[index_df$fw_rv == "fw" & index_df$assay == curr_assay, "primer_num"][1:fw_count[curr_assay][[1]]])
            rv_primers <- c(index_df[index_df$fw_rv == "rv" & index_df$assay == curr_assay, "primer_num"][1:rv_count[curr_assay][[1]]])
          
            if (length(fw_primers) < fw_count[curr_assay][[1]]) {
              stop(
                paste0(
                  "You need at least ",
                  fw_count[curr_assay][[1]],
                  " fw primers. Only found ",
                  length(fw_primers)
                )
              )
            }
            if (length(rv_primers) < rv_count[curr_assay][[1]]) {
              stop(
                paste0(
                  "You need at least ",
                  rv_count[curr_assay][[1]],
                  " rv primers. Only found ",
                  length(rv_primers)
                )
              )
            }
            
            # `plate_count_per_fws` and `curr_fw_plate_count`
            # are used to track the number of plates per fw primer
            # We are tracking this because we increment the rv primer
            # before we start incrementing the fw primer
            plate_count_per_fws <- 1
            curr_fw_plate_count <- 1
    
            # These next four variables are used to splice
            # the fw and rv primer vectors
            fw_left  <- 1
            fw_right <- plate_width
            rv_left  <- 1
            rv_right <- plate_height
    
            sample_i <- 0
            stop     <- FALSE
            for (plate_num in (1:plate_count[curr_assay][[1]])) {
                # We refresh the `well_col` for a new plate
                # It increases for each fw primer
                well_col <- 0
      
                for (fw in fw_primers[fw_left:fw_right]) {
                    well_row_index <- 0
                    well_col       <- well_col + 1
    
                    for (rv in rv_primers[rv_left:rv_right]) {
                        
                        sample_i <- sample_i + 1
    
                        # Stop if we run out of samples
                        if (sample_i > sample_count[curr_assay][[1]]) {
                            stop <- TRUE
                            break
                        }
    
                        # The well will be a string containing a
                        # character and a number; E.g., "A1"
                        well_row_index <- well_row_index + 1
                        well           <- paste0(
                            well_rows[well_row_index],
                            well_col
                        )
    
                        # Add info to metadata
                        meta_row <- which(
                            meta_df$samp_name == sample_ids[curr_assay][[1]][sample_i] &
                            meta_df$assay == curr_assay
                        )[1]
                        
                        meta_df[meta_row, "fw_no"]       <- fw
                        meta_df[meta_row, "rv_no"]       <- rv
                        meta_df[meta_row, "fw_tag"]      <- subset(
                            index_df,
                            primer_num == fw & fw_rv == "fw" & assay == curr_assay,
                            select = tags
                        )$tags
                        meta_df[meta_row, "rv_tag"]      <- subset(
                            index_df,
                            primer_num == rv & fw_rv == "rv" & assay == curr_assay,
                            select = tags
                        )$tags
                        meta_df[meta_row, "plate"]       <- plate_num
                        meta_df[meta_row, "well"]        <- well
                        
                        #if (fw == "52F" & rv == "63R") {print(meta_df[meta_row, ] )}
                        
                    }
                    # Stop looping through fw primers. We ran out of samples.
                    if (stop) {
                        break
                    }
                }
    
                # Stop looping through the plates. We ran out of samples.
                if (stop) {
                    break
                }
    
                # We still have rv primers left
                if (rv_right < rv_count[curr_assay][[1]]) {
                    rv_left             <- rv_left + plate_height
                    rv_right            <- rv_right + plate_height
                    
                    curr_fw_plate_count <- curr_fw_plate_count + 1
    
                # We ran out of rv primers
                } else {
                    # We aren't up to the next fw primer yet
                    if (curr_fw_plate_count < plate_count_per_fws) {
                        
                        curr_fw_plate_count <- curr_fw_plate_count + 1
    
                    # We are up to the next fw primer
                    } else {
                        
                        curr_fw_plate_count <- 1
                        fw_left  <- fw_left + plate_width
                        fw_right <- fw_right + plate_width
                        rv_left  <- 1
                        rv_right <- plate_height
                    }
                }
            }
        }
    } else if (strategy == "UDI") {
        for (curr_assay in assays) {
            # Get vectors of the primers
            num_of_sams <- nrow(meta_df[meta_df$assay == curr_assay, ])
            if (length(index_df[index_df$fw_rv == "FW" & index_df$assay == curr_assay, "primer_num"]) < num_of_sams) {
                stop(
                    paste0(
                        "You need at least ",
                        num_of_sams,
                        " fw primers. Only found ",
                        length(index_df[index_df$fw_rv == "fw" & index_df$assay == curr_assay, "primer_num"])
                    )
                )
            }
            if (length(index_df[index_df$fw_rv == "rv" & index_df$assay == curr_assay, "primer_num"]) < num_of_sams) {
                stop(
                    paste0(
                        "You need at least ",
                        num_of_sams,
                        " rv primers. Only found ",
                        length(index_df[index_df$fw_rv == "rv" & index_df$assay == curr_assay, "primer_num"])
                    )
                )
            }
            fw_primers <- c(index_df[index_df$fw_rv == "fw" & index_df$assay == curr_assay, "primer_num"][1:num_of_sams])
            rv_primers <- c(index_df[index_df$fw_rv == "rv" & index_df$assay == curr_assay, "primer_num"][1:num_of_sams])
            
            
            well_row_index    <- 0
            well_col          <- 1
            plate_num         <- 1
            samp_num_in_plate <- 0
            samp_num_in_col   <- 0
            for (i in (1:length(sample_ids[curr_assay][[1]]))) {
                samp_num_in_plate <- samp_num_in_plate + 1
                samp_num_in_col   <- samp_num_in_col + 1
                # new col
                if (samp_num_in_col > plate_height) {
                    samp_num_in_col <- 1
                    well_row_index <- 0
                    well_col       <- well_col + 1
                }
                # new plate
                if (samp_num_in_plate > sample_count_per_plate + control_count) {
                    samp_num_in_plate <- 1
                    plate_num         <- plate_num + 1
                    well_row_index <- 0
                    well_col          <- 1
                }
                
                sam <- sample_ids[curr_assay][[1]][i]
                fw  <- fw_primers[i]
                rv  <- rv_primers[i]
                
                # Stop if we run out of samples
                if (i > sample_count[curr_assay][[1]]) {
                    stop <- TRUE
                    break
                }
                
                # The well will be a string containing a
                # character and a number; E.g., "A1"
                well_row_index <- well_row_index + 1
                well           <- paste0(
                    well_rows[well_row_index],
                    well_col
                )
                
                # Add info to metadata
                meta_row <- which(
                    meta_df$samp_name == sam &
                        meta_df$assay == curr_assay
                )[1]
                
                meta_df[meta_row, "fw_no"]       <- fw
                meta_df[meta_row, "rv_no"]       <- rv
                meta_df[meta_row, "fw_tag"]      <- subset(
                    index_df,
                    primer_num == fw & fw_rv == "fw" & assay == curr_assay,
                    select = tags
                )$tags
                meta_df[meta_row, "rv_tag"]      <- subset(
                    index_df,
                    primer_num == rv & fw_rv == "rv" & assay == curr_assay,
                    select = tags
                )$tags
                meta_df[meta_row, "plate"]       <- plate_num
                meta_df[meta_row, "well"]        <- well
            }
        }
    }

    return(meta_df)
}
