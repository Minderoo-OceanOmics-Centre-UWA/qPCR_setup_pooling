get_possible_counts <- function(sample_ids,
                                controls,
                                plate_height,
                                plate_width) {

    # Okay, so there's lots of variables here.
    # These variables track possible values based on the number of plates
    # The output of the function is:
    # - a vector of sample cutoffs based on based on number of samples per plate
    #   E.g., c(96, 192, 288, 384)
    # - a vector of possible plate counts
    #   E.g., c(1, 2, 3, 4)
    # - a vector of possible fw counts
    #   E.g., c(12, 12, 12, 24)
    # - a vector of possible rv counts
    #   E.g., c(8, 16, 24, 24)
    sample_count_per_plate <- (plate_height * plate_width) - length(controls)
    sample_count           <- length(sample_ids)
    curr_sample_cutoff     <- 0
    plate_count            <- 0
    sample_cutoffs         <- c()
    possible_plate_counts  <- c()
    possible_fw_counts     <- c()
    possible_rv_counts     <- c()
    #
    # TODO: Maybe add a check to make sure these are multiples of plate_width/plate_height
    #       can maybe use stopifnot() https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/stopifnot
    #
    max_fw_count           <- 36
    max_rv_count           <- 24
    curr_fw_count          <- plate_width
    curr_rv_count          <- 0
    plate_count_per_fws    <- 0
    curr_fw_plate_count    <- 1

    # While we still have more samples
    while (curr_sample_cutoff < sample_count) {
        plate_count           <- plate_count + 1
        possible_plate_counts <- append(possible_plate_counts, plate_count)
        curr_sample_cutoff    <- curr_sample_cutoff + sample_count_per_plate
        sample_cutoffs        <- append(sample_cutoffs, curr_sample_cutoff)

        # We still have rv primers left
        if (curr_rv_count < max_rv_count) {
            curr_rv_count       <- curr_rv_count + plate_height
            plate_count_per_fws <- plate_count_per_fws + 1
            curr_fw_plate_count <- curr_fw_plate_count + 1

        # We ran out of rv primers
        } else {
            # We aren't up to the next fw primer yet
            if (curr_fw_plate_count <= plate_count_per_fws) {
                curr_fw_plate_count <- curr_fw_plate_count + 1

            # We are up to the next fw primer
            } else {
                curr_fw_plate_count <- 2
                curr_fw_count       <- curr_fw_count + plate_width
            }
        }
        
        possible_fw_counts <- append(possible_fw_counts, curr_fw_count)
        possible_rv_counts <- append(possible_rv_counts, curr_rv_count)
    }

    list(
        plates = possible_plate_counts,
        fws = possible_fw_counts,
        rvs = possible_rv_counts,
        cutoffs = sample_cutoffs
    )
}
