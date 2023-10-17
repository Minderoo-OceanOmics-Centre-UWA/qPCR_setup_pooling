get_possible_counts <- function(sample_ids, controls, plate_height, plate_width) {
    sample_count_per_plate <- (plate_height * plate_width) - length(controls)
    sample_count <- length(sample_ids)
    curr_sample_cutoff <- 0
    plate_count <- 0
    sample_cutoffs <- c()
    possible_plate_counts <- c()
    possible_fw_counts <- c()
    possible_rv_counts <- c()
    max_fw_count <- 36
    max_rv_count <- 24 # add a check to make sure these are multiples of plate_width/plate_height
    curr_fw_count <- plate_width
    curr_rv_count <- 0
    plate_count_per_fws <- 0
    curr_fw_plate_count <- 1
    while (curr_sample_cutoff < sample_count) {
        plate_count           <- plate_count + 1
        possible_plate_counts <- append(possible_plate_counts, plate_count)
        curr_sample_cutoff       <- curr_sample_cutoff + sample_count_per_plate
        sample_cutoffs           <- append(sample_cutoffs, curr_sample_cutoff)

        if (curr_rv_count < max_rv_count) {
            curr_rv_count <- curr_rv_count + plate_height
            plate_count_per_fws <- plate_count_per_fws + 1
        } else {
            if (curr_fw_plate_count <= plate_count_per_fws) {
                curr_fw_plate_count <- curr_fw_plate_count + 1
            } else {
                curr_fw_plate_count <- 0
                curr_fw_count <- curr_fw_count + plate_width
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