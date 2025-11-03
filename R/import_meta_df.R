import_meta_df <- function(excel_file, run) {
    excel_sheets <- excel_sheets(excel_file)
    meta_sheet   <- grep(
        glob2rx("sampleMetadata"),
        excel_sheets,
        ignore.case = TRUE,
        value = TRUE
    )

    # Validate number of sheets called 'metadata'
    if (length(meta_sheet) == 0) {
        stop("No sheets found nemed 'sampleMetadata'")
    }
    if (length(meta_sheet) > 1) {
        stop("Multiple sheets found named 'sampleMetadata'")
    }

    # Import metadata
    meta_df   <- read_excel(excel_file, sheet = "sampleMetadata") %>%
        as.data.frame()
    
    num_rows_to_skip = 0
    skip_rows = FALSE
    for (row in rownames(meta_df)) {
        num_rows_to_skip = num_rows_to_skip + 1
        if (meta_df[row, 1] == "samp_name") {
            skip_rows = TRUE
            break
        }
    }
    
    if (skip_rows) {
        meta_df   <- read_excel(excel_file, sheet = "sampleMetadata", skip = num_rows_to_skip) %>%
            as.data.frame()
    }

    # Validata metadata columns
    if (!("samp_name" %in% colnames(meta_df))) {
        stop(paste0("Sheet 'sampleMetadata' must contain a 'samp_name' column"))
    }

    # Subset data if there is a sequencing run column
    if ("sequencingrun" %in% colnames(meta_df)) {
        meta_df <- subset(meta_df, sequencingrun == run)
    }

    # Add empty columns
    meta_df["fw_no"]       <- NULL
    meta_df["rv_no"]       <- NULL
    meta_df["fw_tag"]      <- NULL
    meta_df["rv_tag"]      <- NULL
    meta_df["fw_full_seq"] <- NULL
    meta_df["rv_full_seq"] <- NULL
    meta_df["assay"]       <- NULL
    meta_df["plate"]       <- NULL
    meta_df["well"]        <- NULL

    # Get sample ids
    sample_ids             <- meta_df$samp_name

    # Make sure there is no duplicate samples
    if (length(sample_ids) != length(unique(sample_ids))) {
        duplicates <- meta_df[
            duplicated(meta_df$samp_name) |
            duplicated(meta_df$samp_name, fromLast = TRUE),
            "samp_name"
        ]
        stop(
            paste0("You can't have duplicate sample ids. Duplicates: ",
            duplicates)
        )
    }

    return(meta_df)
}
