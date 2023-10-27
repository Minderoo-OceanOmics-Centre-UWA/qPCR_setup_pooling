print_failed_rep_message <- function(rep_failed_summary, assays) {
    cat("\n")
    cat("--------------------------------------------------------------------------------")
    cat("\nFailed Replicate Message: ")
    for (curr_assay in assays) {
        count_eq_one <- rep_failed_summary %>%
            filter(assay == curr_assay, count_discard == 1) %>%
            n_distinct()
        count_ge_two <- rep_failed_summary %>%
            filter(assay == curr_assay, count_discard >= 2) %>%
            n_distinct()
        cat(
            (ifelse(count_eq_one == 1, "\nThere is", "\nThere are")),
            count_eq_one,
            (ifelse(count_eq_one == 1, "sample in", "samples in")),
            curr_assay,
            "assay with 1 replicate to be discarded",
            (ifelse(count_ge_two == 1, "\nThere is", "\nThere are")),
            count_ge_two,
            (ifelse(count_ge_two == 1, "sample in", "samples in")),
            curr_assay,
            "assay with 2 or more failed replicates\n"
        )
    }
    #cat("\n")
    cat("--------------------------------------------------------------------------------")
    cat("\n\n")
}