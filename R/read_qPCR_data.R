# Loop to manipulate all files simultaneously
read_qPCR_data <- function(file, assays, plate_numbers) {
    # loads in data
    data <- read.delim(file)

    # get info from sample source, assay and plate id details from file name
    data$fileName <- gsub(
        ".txt",
        "",
        str_split(file, "\\/\\/", simplify = TRUE)[, 2]
    )
    description       <- strsplit(data$fileName, "_")
    desc_count        <- length(description[[1]])
    data$assay        <- sapply(description, "[", (desc_count - 1))
    data$plate_number <- sapply(description, "[", desc_count)

    curr_assay <- unique(data$assay)
    curr_plate <- unique(data$plate_number)
    if (! curr_assay %in% assays) {
        stop(
            paste0(
                "Error: assay - ",
                curr_assay,
                " taken from file ",
                file,
                " not found in ",
                assays,
                " Make sure filename ends in _$assay_$plate.txt"
            )
        )
    }

    if (! curr_plate %in% plate_numbers) {
        stop(
            paste0(
                "Error: plate - ",
                curr_plate,
                " taken from file ",
                file,
                " not found in ",
                plate_numbers,
                " Make sure filename ends in _$assay_$plate.txt"
            )
        )
    }

    return(data)
}
