import_index_df <- function(excel_file, assays) {
    excel_sheets <- excel_sheets(excel_file)

    index_df     <- data.frame(
        primer_num = character(),
        tags = character(),
        fw_rv = character(),
        stringsAsFactors = FALSE
    )

    for (assay in assays) {
        # Get sheet names that match the correct pattern
        curr_index_sheet <- grep(
            glob2rx(
                paste0("*", assay, "_index*|*index_", assay, "*")
            ),
            excel_sheets,
            ignore.case = TRUE,
            value = TRUE
        )

        # Validate number of sheets that matched pattern
        if (length(curr_index_sheet) == 0) {
            stop(
                paste0(
                    "No sheets found matching pattern: ",
                    "*", assay, "_index*|*index_", assay, "*"
                )
            )
        }
        if (length(curr_index_sheet) > 1) {
            stop(
                paste0(
                    "Multiple sheets found matching pattern: ",
                    "*", assay, "_index*|*index_", assay, "*"
                )
            )
        }

        # Import index sheet
        curr_index_df <- read_excel(excel_file, sheet = curr_index_sheet) %>%
            as.data.frame()

        curr_index_df <- curr_index_df %>% dplyr::select(-contains("..."))
        
        # Validate index sheet columns
        if (!("primer_num" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'primer_num' column"
                )
            )
        }
        if (!("tags" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'tags' column"
                )
            )
        }
        if (!("fw_rv" %in% colnames(curr_index_df))) {
            stop(
                paste0(
                    "Sheet '",
                    curr_index_sheet,
                    "' must contain a 'fw_rv' column"
                )
            )
        }
        if (length(curr_index_df$tags) !=
        length(unique(curr_index_df$tags))) {
            duplicates <- curr_index_df[
                duplicated(curr_index_df$tags) |
                duplicated(curr_index_df$tags, fromLast = TRUE),
                "tags"
            ]
            stop(
                paste0(
                    "You can't have duplicate tags. Duplicates: ",
                    duplicates
                )
            )
        }

        # Make sure fw_rv column only has 'fw' or 'rv' values
        unique_fw_rv       <- unique(curr_index_df$fw_rv)
        for (i in unique_fw_rv) {
            if (i != "fw" && i != "rv") {
                stop(
                    paste0(
                        "You can't have ",
                        i,
                        " in fw_rv column. Valid values are 'fw' or 'rv'"
                    )
                )
            }
        }

        curr_index_df$assay      <- assay
        curr_index_df$primer_seq <- NULL
        index_df                 <- rbind(index_df, curr_index_df)
    }

    return(index_df)
}
