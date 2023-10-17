reformat_meta_df <- function(meta_df,
                             index_df,
                             assays,
                             controls,
                             plate_count,
                             fw_count,
                             rv_count,
                             plate_height,
                             plate_width) {

    control_count <- length(controls)
    sample_count_per_plate <- (plate_height * plate_width) - control_count

    # add sample_controls
    left_of_controls <- sample_count_per_plate
    for (plate_num in 1:plate_count) {

        control_rows <- data.frame(
            SAMPLEID = paste0(controls, "_", plate_num),
            SEQUENCINGRUN = run
        )

        if (left_of_controls > nrow(meta_df)) {
            left_of_controls <- nrow(meta_df)
            meta_df <- rbind(meta_df, control_rows)
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

    # Get a vector of just the primers we need
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

    # Track the well we are up to so we can include the well info
    well_rows      <- c()
    char           <- "A"
    for (i in 1:plate_height) {
        well_rows <- append(well_rows, char)
        # Get next letter
        char      <- LETTERS[which(LETTERS == char) + 1L]
    }

    # There is 12 columns and 8 rows per plate,
    # so these for loops will iterate 96 times
    for (assay in assays) {
        sample_i <- 0
        plate_count_per_fws <- 0
        curr_fw_plate_count <- 1
        fw_left  <- 1
        fw_right <- plate_width
        rv_left  <- 1
        rv_right <- plate_height
        stop <- FALSE
        for (plate_num in (1:plate_count)) {
            well_col <- 0
            for (fw in fw_primers[fw_left:fw_right]) {
                well_row_index <- 0
                well_col       <- well_col + 1

                for (rv in rv_primers[rv_left:rv_right]) {
                    sample_i <- sample_i + 1
                    if (sample_i > sample_count) {
                        stop <- TRUE
                        break
                    }
                    well_row_index <- well_row_index + 1
                    well           <- paste0(
                        well_rows[well_row_index],
                        well_col
                    )

                    # Add primer info to metadata
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
                if (stop) {
                    break
                }
            }
            if (stop) {
                break
            }

            if (rv_right < rv_count) {
                rv_left <- rv_left + plate_height
                rv_right <- rv_right + plate_height
                plate_count_per_fws <- plate_count_per_fws + 1
            } else {
                if (curr_fw_plate_count <= plate_count_per_fws) {
                    curr_fw_plate_count <- curr_fw_plate_count + 1
                } else {
                    curr_fw_plate_count <- 0
                    fw_left <- fw_left + plate_width
                    fw_right <- fw_right + plate_width
                }
            }
        }
    }

    return(meta_df)
}
