reformat_meta_df <- function(meta_df,
                             index_df,
                             assays,
                             controls,
                             plate_count,
                             fw_count,
                             rv_count,
                             plate_height,
                             plate_width,
                             run) {
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
                             sample_ids) {

    control_count <- length(controls)
    sample_count_per_plate <- (plate_height * plate_width) - control_count
    sample_rows <- list()
    sample_count <- list()
    
    for (assay in assays) {
        sample_rows[assay][[1]] <- data.frame(
            SAMPLEID = sample_ids[assay][[1]],
            SEQUENCINGRUN = run
        )
        
        # Add control samples to meta_df.
        # We need to add the control samples at the end of each plate.
        left_of_controls <- sample_count_per_plate
        fake_plates <- 0
        for (plate_num in 1:plate_count[assay][[1]]) {
          if (plate_num == 1) {
              if (all(grepl("^FAKESAMPLE_", sample_rows[assay][[1]][1:left_of_controls, ]$SAMPLEID))) {
                  fake_plates <- fake_plates + 1
              }  
          } else {
            if (all(grepl("^FAKESAMPLE_", sample_rows[assay][[1]][(prev_left_of_controls + control_count + 1):left_of_controls, ]$SAMPLEID))) {
              fake_plates <- fake_plates + 1
            }  
          }
          
          if ((plate_num - fake_plates) == 0) {
            control_rows <- data.frame(
              SAMPLEID = paste0("FAKESAMPLE_", controls, "_", plate_num),
              SEQUENCINGRUN = run
            )
          } else {
            control_rows <- data.frame(
              SAMPLEID = paste0(controls, "_", (plate_num - fake_plates)),
              SEQUENCINGRUN = run
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
        
        sample_ids[assay][[1]] <- sample_rows[assay][[1]]$SAMPLEID
        sample_count[assay][[1]] <- length(sample_ids[assay][[1]])
    }

    

    # We need to duplicate the samples for each assay
    sample_rows <- rbind.fill(sample_rows)
    assay_vec   <- c()
    for (assay in assays) {
        assay_vec <- c(assay_vec, rep(assay, length.out = sample_count[assay][[1]]))
    }

    sample_rows$ASSAY      <- assay_vec
    row.names(sample_rows) <- NULL
    
    # Create a vector of characters based on plate height.
    well_rows      <- c()
    char           <- "A"
    for (i in 1:plate_height) {
      well_rows <- append(well_rows, char)
      # Get next letter
      char      <- LETTERS[which(LETTERS == char) + 1L]
    }
    
    meta_df <- merge(sample_rows, meta_df, by = c("SAMPLEID", "SEQUENCINGRUN"), all.x = TRUE)

    meta_df <- meta_df[match(paste(sample_rows$SAMPLEID, sample_rows$ASSAY), 
                                    paste(meta_df$SAMPLEID, meta_df$ASSAY)), ]
    rownames(meta_df) <- NULL
    
    
    #tmp_meta_df <- meta_df
    
    
    #meta_df <- tmp_meta_df
    
    for (assay in assays) {
        # Get vectors of the primers
        fw_primers <- c(index_df[index_df$FWRV == "FW" & index_df$ASSAY == assay, "PRIMERNUM"][1:fw_count[assay][[1]]])
        rv_primers <- c(index_df[index_df$FWRV == "RV" & index_df$ASSAY == assay, "PRIMERNUM"][1:rv_count[assay][[1]]])
      
        if (length(fw_primers) < fw_count[assay][[1]]) {
          stop(
            paste0(
              "You need at least ",
              fw_count[assay][[1]],
              " fw primers. Only found ",
              length(fw_primers)
            )
          )
        }
        if (length(rv_primers) < rv_count[assay][[1]]) {
          stop(
            paste0(
              "You need at least ",
              rv_count[assay][[1]],
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
        for (plate_num in (1:plate_count[assay][[1]])) {
            # We refresh the `well_col` for a new plate
            # It increases for each fw primer
            well_col <- 0
  
            for (fw in fw_primers[fw_left:fw_right]) {
                well_row_index <- 0
                well_col       <- well_col + 1
                
                
                #print(meta_df[meta_df$SAMPLEID == "ABv9_DC_3", ])

                for (rv in rv_primers[rv_left:rv_right]) {
                    
                    sample_i <- sample_i + 1

                    # Stop if we run out of samples
                    if (sample_i > sample_count[assay][[1]]) {
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
                        meta_df$SAMPLEID == sample_ids[assay][[1]][sample_i] &
                        meta_df$ASSAY == assay
                    )[1]
                    
                    meta_df[meta_row, "FW_NO"]       <- fw
                    meta_df[meta_row, "RV_NO"]       <- rv
                    meta_df[meta_row, "FW_TAG"]      <- subset(
                        index_df,
                        PRIMERNUM == fw & FWRV == "FW" & ASSAY == assay,
                        select = TAGS
                    )$TAGS
                    meta_df[meta_row, "RV_TAG"]      <- subset(
                        index_df,
                        PRIMERNUM == rv & FWRV == "RV" & ASSAY == assay,
                        select = TAGS
                    )$TAGS
                    meta_df[meta_row, "PLATE"]       <- plate_num
                    meta_df[meta_row, "WELL"]        <- well
                    
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
            if (rv_right < rv_count[assay][[1]]) {
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

    return(meta_df)
}

    control_count <- length(controls)
    sample_count_per_plate <- (plate_height * plate_width) - control_count

    # Add control samples to meta_df.
    # We need to add the control samples at the end of each plate.
    left_of_controls <- sample_count_per_plate
    for (plate_num in 1:plate_count) {
        control_rows <- data.frame(
            SAMPLEID = paste0(controls, "_", plate_num),
            SEQUENCINGRUN = run
        )

        # This will be true for the last plate when there isn't enough
        # samples to fill the plate
        if (left_of_controls > nrow(meta_df)) {
            left_of_controls <- nrow(meta_df)
            meta_df <- rbind.fill(meta_df, control_rows)

        # This will be true for most plates
        } else {
            meta_df <- rbind.fill(
                meta_df[1:left_of_controls, ],
                control_rows,
                meta_df[(left_of_controls + 1):nrow(meta_df), ]
            )
        }

        left_of_controls <- (
            left_of_controls +
            control_count +
            sample_count_per_plate
        )
    }

    sample_ids <- meta_df$SAMPLEID
    sample_count <- length(sample_ids)

    # We need to duplicate the samples in the metadata for each assay
    meta_df <- meta_df[
        rep(
            seq_len(
                nrow(meta_df)),
                each = length(assays)
            ),
    ]

    meta_df$ASSAY      <- rep(assays, length.out = nrow(meta_df))
    row.names(meta_df) <- NULL

    # Get vectors of the primers
    fw_primers <- c(index_df[index_df$FWRV == "FW" & index_df$ASSAY == assays[1], "PRIMERNUM"][1:fw_count])
    rv_primers <- c(index_df[index_df$FWRV == "RV" & index_df$ASSAY == assays[1], "PRIMERNUM"][1:rv_count])

    if (length(fw_primers) < fw_count) {
        stop(
            paste0(
                "You need at least ",
                fw_count,
                " fw primers. Only found ",
                length(fw_primers)
            )
        )
    }
    if (length(rv_primers) < rv_count) {
        stop(
            paste0(
                "You need at least ",
                rv_count,
                " rv primers. Only found ",
                length(rv_primers)
            )
        )
    }

    # Create a vector of characters based on plate height.
    well_rows      <- c()
    char           <- "A"
    for (i in 1:plate_height) {
        well_rows <- append(well_rows, char)
        # Get next letter
        char      <- LETTERS[which(LETTERS == char) + 1L]
    }

    for (assay in assays) {
        # Get vectors of the primers
        fw_primers <- c(index_df[index_df$FWRV == "FW" & index_df$ASSAY == assay, "PRIMERNUM"][1:fw_count])
        rv_primers <- c(index_df[index_df$FWRV == "RV" & index_df$ASSAY == assay, "PRIMERNUM"][1:rv_count])
      
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
        for (plate_num in (1:plate_count)) {
            # We refresh the `well_col` for a new plate
            # It increases for each fw primer
            well_col <- 0
  
            for (fw in fw_primers[fw_left:fw_right]) {
                well_row_index <- 0
                well_col       <- well_col + 1
                
                
                #print(meta_df[meta_df$SAMPLEID == "ABv9_DC_3", ])

                for (rv in rv_primers[rv_left:rv_right]) {
                    sample_i <- sample_i + 1

                    # Stop if we run out of samples
                    if (sample_i > sample_count) {
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
                        meta_df$SAMPLEID == sample_ids[sample_i] &
                        meta_df$ASSAY == assay
                    )[1]
                    meta_df[meta_row, "FW_NO"]       <- fw
                    meta_df[meta_row, "RV_NO"]       <- rv
                    meta_df[meta_row, "FW_TAG"]      <- subset(
                        index_df,
                        PRIMERNUM == fw & FWRV == "FW" & ASSAY == assay,
                        select = TAGS
                    )$TAGS
                    meta_df[meta_row, "RV_TAG"]      <- subset(
                        index_df,
                        PRIMERNUM == rv & FWRV == "RV" & ASSAY == assay,
                        select = TAGS
                    )$TAGS
                    meta_df[meta_row, "PLATE"]       <- plate_num
                    meta_df[meta_row, "WELL"]        <- well
                    
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
            if (rv_right < rv_count) {
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

    return(meta_df)
}
