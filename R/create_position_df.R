create_position_df <- function(meta_df,
                               assays,
                               plate_count,
                               plate_height,
                               plate_width,
                               control_pattern) {

    sample_ids <- unique(meta_df$SAMPLEID)
    sample_count <- length(sample_ids)

    # Create a position_df to hold information that can link
    # the samples to the big plate positions
    position_df                     <- meta_df[, c("SAMPLEID", "ASSAY")]
    colnames(position_df)[colnames(position_df) == "SAMPLEID"] <- "sample_id"
    position_df["Pos"]              <- NULL
    position_df["sample_replicate"] <- NULL
    position_df["replicate"]        <- NULL
    position_df["plate_number"]     <- NULL
    position_df["assay_plate_number_Pos"] <- NULL
    position_df["plate_number_Pos"] <- NULL
    position_df["sample_type"]      <- NULL
    replicates                      <- c("1", "2", "3", "pool")
    position_df <- position_df[
        rep(
            seq_len(
                nrow(
                    position_df
                )
            ),
            each = length(replicates)
        ),
    ]
    position_df$replicate <- rep(
        replicates,
        length.out = nrow(position_df)
    )
    row.names(position_df)          <- NULL

    well_rows      <- c()
    char           <- "A"
    for (i in 1:(plate_height*2)) {
        well_rows <- append(well_rows, char)
        # Get next letter
        char      <- LETTERS[which(LETTERS == char) + 1L]
    }

    for (assay in assays) {
        sample_i <- 1
        for (plate_num in (1:plate_count)) {
            # Track the well and replicate we are up to \
            well_row_index <- 1
            well_col       <- 1
            replicate      <- "1"

            # There is 24 columns and 16 rows per big plate, 
            # so we need to loop 384 times
            for (i in 1:((plate_height*2)*(plate_width*2))) {
                well <- paste0(well_rows[well_row_index], well_col)


                if (sample_i > sample_count) {
                    break
                }
                curr_sample     <- sample_ids[sample_i]

                # Samples with these patterns should be flaged as control samples
                if (grepl(control_pattern, c(curr_sample))) {
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
                    well_col <- well_col + 1

                # If the current replicate is '2', we need to move one column to the left and one row down
                } else if (replicate == "2") {
                    replicate      <- "3"
                    well_row_index <- well_row_index + 1
                    well_col <- well_col - 1

                # If the current replicate is '3', we need to move one column to the right
                } else if (replicate == "3") {
                    replicate      <- "pool"
                    well_col <- well_col + 1
                        
                    # If the current replicate is 'pool', we need to check if we are up to the last row 
                } else if (replicate == "pool") {
                replicate      <- "1"
                sample_i <- sample_i + 1 

                    # If it is the last row, move to the first row, and move across one columns
                    if (well_row_index == 16) {
                        well_row_index <- 1
                        well_col <- well_col + 1

                    # If it is not the last row, move down one row, but stay in the same column    
                    } else {
                        well_row_index <- well_row_index + 1
                        well_col <- well_col - 1
                            
                    }
                }
            }
        }
    }

    return(position_df)
}