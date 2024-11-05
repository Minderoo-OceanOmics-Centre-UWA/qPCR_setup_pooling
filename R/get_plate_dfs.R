get_plate_dfs <- function(first_fw,
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
                          meta_df) {

    sample_ids <- meta_df$SAMPLEID[meta_df$ASSAY == assay]
    sample_count <- length(sample_ids)

    # Get vectors of the primers
    fw_primers <- c(index_df[index_df$FWRV == "FW" & index_df$ASSAY == assay, "PRIMERNUM"][1:fw_count])
    rv_primers <- c(index_df[index_df$FWRV == "RV" & index_df$ASSAY == assay, "PRIMERNUM"][1:rv_count])

    plate     <- data.frame(
        matrix(
            ncol = plate_width,
            nrow = plate_height
        )
    )
    big_plate <- data.frame(
        matrix(
            ncol = (plate_width * 2),
            nrow = (plate_height * 2)
        )
    )

    # We need to duplicate the primers for the big_plates,
    # but each value must be unique, so let's append _1, or _2
    doubled_fw_primers <- rep(fw_primers[first_fw:last_fw], each = 2)
    doubled_fw_primers <- paste0(doubled_fw_primers, c("_1", "_2"))
    doubled_rv_primers <- rep(rv_primers[first_rv:last_rv], each = 2)
    doubled_rv_primers <- paste0(doubled_rv_primers, c("_1", "_2"))

    colnames(plate)     <- fw_primers[first_fw:last_fw]
    colnames(big_plate) <- doubled_fw_primers
    rownames(plate)     <- rv_primers[first_rv:last_rv]
    rownames(big_plate) <- doubled_rv_primers

    stop <- FALSE

    # Loop through each cell of small plate
    sample_i <- first_sample
    for (col in colnames(plate)) {
        for (row in rownames(plate)) {

            # If we run out of samples
            if (sample_i > sample_count) {
                stop <- TRUE
                break
            }
            curr_sample     <- sample_ids[sample_i]
            
            if (grepl("^FAKESAMPLE", curr_sample)) {
              plate[row, col] <- NA
              big_plate[paste0(row, "_1"), paste0(col, "_1")] <- NA
              big_plate[paste0(row, "_1"), paste0(col, "_2")] <- NA
              big_plate[paste0(row, "_2"), paste0(col, "_1")] <- NA
              big_plate[paste0(row, "_2"), paste0(col, "_2")] <- NA
            } else {
              # Add the current sample to the current cell
              plate[row, col] <- curr_sample
              
              # Each sample should be added to the big plates four times
              big_plate[paste0(row, "_1"), paste0(col, "_1")] <- paste0(
                curr_sample, "-1"
              )
              big_plate[paste0(row, "_1"), paste0(col, "_2")] <- paste0(
                curr_sample, "-2"
              )
              big_plate[paste0(row, "_2"), paste0(col, "_1")] <- paste0(
                curr_sample, "-3"
              )
              big_plate[paste0(row, "_2"), paste0(col, "_2")] <- paste0(
                curr_sample, "-pool"
              ) 
            }

            sample_i <- sample_i + 1
        }

        # Stop. We ran out of samples.
        if (stop) {
            break
        }
    }

    list(small_plate = plate, big_plate = big_plate)
}
