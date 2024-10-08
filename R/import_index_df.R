import_index_df <- function(excel_file, assays) {
    excel_sheets <- excel_sheets(excel_file)

    index_df     <- data.frame(
        PRIMERNUM = character(),
        TAGS = character(),
        FWRV = character(),
        stringsAsFactors = FALSE
    )

    for (assay in assays) {
        # Get sheet names that match the correct pattern
        curr_index_sheet <- grep(
            glob2rx(
                paste0("*", assay, "*index*|*index*", assay, "*")
            ),
            excel_sheets,
            ignore.case = TRUE,
            value = TRUE
        )

        # Validate number of sheets that matched pattern
        if (length(curr_index_sheet) == 0) {
            stop(
                paste0(
                    "No sheets found containing words 'index' and ",
                    assay
                )
            )
        }
        if (length(curr_index_sheet) > 1) {
            stop(
                paste0(
                    "Multiple sheets found containing words 'index' and ",
                    assay
                )
            )
        }

        # Import index sheet
        curr_index_df <- read_excel(excel_file, sheet = curr_index_sheet) %>%
            as.data.frame()
        names(curr_index_df) <- gsub(" ", "", names(curr_index_df)) %>% {
                gsub("_", "", .)
            } %>%
            toupper()

        curr_index_df <- curr_index_df %>% dplyr::select(-contains("..."))
        
        # Validate index sheet columns
        if (!("PRIMERNUM" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'primer_num' column"
                )
            )
        }
        if (!("TAGS" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'tags' column"
                )
            )
        }
        if (!("FWRV" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'fw_rv' column"
                )
            )
        }
        if (length(curr_index_df$TAGS) !=
        length(unique(curr_index_df$TAGS))) {
            duplicates <- curr_index_df[
                duplicated(curr_index_df$TAGS) |
                duplicated(curr_index_df$TAGS, fromLast = TRUE),
                "TAGS"
            ]
            stop(
                paste0(
                    "You can't have duplicate tags. Duplicates: ",
                    duplicates
                )
            )
        }

        # Make sure fw_rv column only has 'fw' or 'rv' values
        curr_index_df$FWRV <- toupper(curr_index_df$FWRV)
        unique_fw_rv       <- unique(curr_index_df$FWRV)
        for (i in unique_fw_rv) {
            if (i != "FW" && i != "RV") {
                stop(
                    paste0(
                        "You can't have ",
                        i,
                        " in fw_rv column. Valid values are 'fw' or 'rv'"
                    )
                )
            }
        }

        curr_index_df$ASSAY <- assay
        index_df            <- rbind(index_df, curr_index_df)
    }

    return(index_df)
}
