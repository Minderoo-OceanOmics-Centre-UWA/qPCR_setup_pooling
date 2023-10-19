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
            meta_df <- rbind(meta_df, control_rows)

        # This will be true for most plates
        } else {
            meta_df <- rbind(
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
    fw_primers <- c(index_df[index_df$FWRV == "FW", "PRIMERNUM"][1:fw_count])
    rv_primers <- c(index_df[index_df$FWRV == "RV", "PRIMERNUM"][1:rv_count])

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
        # `plate_count_per_fws` and `curr_fw_plate_count`
        # are used to track the number of plates per fw primer
        # We are tracking this because we increment the rv primer
        # before we start incrementing the fw primer
        plate_count_per_fws <- 0
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
                    )
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
                    meta_df[meta_row, "FW_FULL_SEQ"] <- subset(
                        index_df,
                        PRIMERNUM == fw & FWRV == "FW" & ASSAY == assay,
                        select = PRIMERSEQ
                    )$PRIMERSEQ
                    meta_df[meta_row, "RV_FULL_SEQ"] <- subset(
                        index_df,
                        PRIMERNUM == rv & FWRV == "RV" & ASSAY == assay,
                        select = PRIMERSEQ
                    )$PRIMERSEQ
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

            #
            # TODO: I may need to test this more to make sure it works
            #
            # We still have rv primers left
            if (rv_right < rv_count) {
                rv_left             <- rv_left + plate_height
                rv_right            <- rv_right + plate_height
                plate_count_per_fws <- plate_count_per_fws + 1

            # We ran out of rv primers
            } else {
                # We aren't up to the next fw primer yet
                if (curr_fw_plate_count <= plate_count_per_fws) {
                    curr_fw_plate_count <- curr_fw_plate_count + 1

                # We are up to the next fw primer
                } else {
                    curr_fw_plate_count <- 0
                    fw_left  <- fw_left + plate_width
                    fw_right <- fw_right + plate_width
                }
            }
        }
    }

    return(meta_df)
}
