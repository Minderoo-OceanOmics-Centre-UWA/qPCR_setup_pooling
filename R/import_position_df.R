import_position_df <- function(excel_file) {
    excel_sheets   <- excel_sheets(excel_file)
    position_sheet <- grep(
        glob2rx("position_df"),
        excel_sheets,
        ignore.case = TRUE,
        value = TRUE
    )

    # Validate number of sheets called 'position_df'
    if (length(position_sheet) == 0) {
        stop("No sheets found named 'position_df'")
    }
    if (length(position_sheet) > 1) {
        stop("Multiple sheets found named 'position_df'")
    }

    # Import df
    position_df   <- read_excel(excel_file, sheet = "position_df") %>%
        as.data.frame()

    # Validata columns
    if (!("assay" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain an 'assay' column")
    }
    if (!("plate_number" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain a 'plate_number' column")
    }
    if (!("plate_number_pos" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain a 'plate_number_pos' column")
    }
    if (!("replicate" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain a 'replicate' column")
    }
    if (!("sample_replicate" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain a 'sample_replicate' column")
    }
    if (!("sample_type" %in% colnames(position_df))) {
        stop("Sheet 'position_df' must contain a 'sample_type' column")
    }

    return(position_df)
}