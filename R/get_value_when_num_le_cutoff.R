get_value_when_num_le_cutoff <- function(num, cutoffs, values) {
    i <- 1

    for (cutoff in cutoffs) {
        if (num <= cutoff) {
            return(values[i])
        }
        i <- i + 1
    }
    stop(
        paste0(
            num, " is greater than ", cutoffs[-1]
        )
    )
}
