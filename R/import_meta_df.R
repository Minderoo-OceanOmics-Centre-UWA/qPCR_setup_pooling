import_meta_df <- function(excel_file, run) {
    excel_sheets <- excel_sheets(excel_file)
    meta_sheet   <- grep(
        glob2rx("metadata"),
        excel_sheets,
        ignore.case = TRUE,
        value = TRUE
    )

    # Validate number of sheets called 'metadata'
    if (length(meta_sheet) == 0) {
        stop("No sheets found nemed 'metadata'")
    }
    if (length(meta_sheet) > 1) {
        stop("Multiple sheets found named 'metadata'")
    }

    # Import metadata
    meta_df   <- read_excel(excel_file, sheet = "metadata") %>%
        as.data.frame()
    names(meta_df) <- gsub(" ", "", names(meta_df)) %>% {
            gsub("_", "", .)
        } %>%
        toupper()

    # Validata metadata columns
    if (!("SAMPLEID" %in% colnames(meta_df))) {
        stop(paste0("Sheet 'metadata' must contain a 'sample_id' column"))
    }

    # Subset data if there is a sequencing run column
    if ("SEQUENCINGRUN" %in% colnames(meta_df)) {
        meta_df <- subset(meta_df, SEQUENCINGRUN == run)
    }

    # Add empty columns
    meta_df["FW_NO"]       <- NULL
    meta_df["RV_NO"]       <- NULL
    meta_df["FW_TAG"]      <- NULL
    meta_df["RV_TAG"]      <- NULL
    meta_df["FW_FULL_SEQ"] <- NULL
    meta_df["RV_FULL_SEQ"] <- NULL
    meta_df["ASSAY"]       <- NULL
    meta_df["PLATE"]       <- NULL
    meta_df["WELL"]        <- NULL

    # Get sample ids
    sample_ids             <- meta_df$SAMPLEID

    # Make sure there is no duplicate samples
    if (length(sample_ids) != length(unique(sample_ids))) {
        duplicates <- meta_df[
            duplicated(meta_df$SAMPLE_ID) |
            duplicated(meta_df$SAMPLE_ID, fromLast = TRUE),
            "SAMPLE_ID"
        ]
        stop(
            paste0("You can't have duplicate sample ids. Duplicates: ",
            duplicates)
        )
    }

    return(meta_df)
}