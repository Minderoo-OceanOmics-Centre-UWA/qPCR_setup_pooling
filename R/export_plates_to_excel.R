export_plates_to_excel <- function(assays,
                                   plates,
                                   big_plates,
                                   meta_df,
                                   position_df,
                                   output_file,
                                   plate_height) {
    # Create a workbook object
    wb <- createWorkbook()

    # Add a sheet to the workbook
    addWorksheet(wb, "plates")
    addWorksheet(wb, "big_plates")
    addWorksheet(wb, "metadata")
    addWorksheet(wb, "position_df")

    # We will loop through the assays and plate numbers and use those values as
    # keys to get the correct data for populating the output Excel file
    # 'row_start' and 'big_row_start' are to track the rows we are up to in the
    # Excel file so that the next plate is placed under the previous plate in
    # the Excel sheet
    row_start        <- 1
    big_row_start    <- 1
    plate_count      <- length(plates) / length(assays)

    for (assay in assays) {
        for (plate_num in 1:plate_count) {
            key <- paste0(assay, plate_num)
            plates[key][[1]] <- rownames_to_column(
                plates[key][[1]], var = paste0(assay, " plate ", plate_num)
            )

            big_plates[key][[1]] <- rownames_to_column(
                big_plates[key][[1]],
                var = paste0(assay, " plate ", plate_num)
            )

            big_plates[key][[1]][, 1] <- gsub(
                "_[12]$",
                "",
                big_plates[key][[1]][, 1]
            )

            colnames(big_plates[key][[1]])[-1] <- gsub(
                "_[12]$",
                "",
                colnames(big_plates[key][[1]])[-1]
            )

            plates[key][[1]] <- rbind(
                colnames(plates[key][[1]]),
                plates[key][[1]]
            )

            colnames(plates[key][[1]]) <- ""

            big_plates[key][[1]] <- rbind(
                colnames(big_plates[key][[1]]),
                big_plates[key][[1]]
            )

            colnames(big_plates[key][[1]]) <- ""

            plates[key][[1]] <- replace(
                plates[key][[1]],
                is.na(plates[key][[1]]) | plates[key][[1]] == "NA",
                ""
            )

            big_plates[key][[1]] <- replace(
                big_plates[key][[1]],
                is.na(big_plates[key][[1]]) | big_plates[key][[1]] == "NA",
                ""
            )

            # Write the data to the sheet
            writeData(
                wb,
                sheet = "plates",
                x = plates[key][[1]],
                startRow = row_start,
                startCol = 1,
                colNames = FALSE,
                rowNames = FALSE
            )

            writeData(
                wb,
                sheet = "big_plates",
                x = big_plates[key][[1]],
                startRow = big_row_start,
                startCol = 1,
                colNames = FALSE,
                rowNames = FALSE
            )

            row_start     <- row_start + plate_height + 2
            big_row_start <- big_row_start + (plate_height * 2) + 2
        }
    }

    # We can't forget the add the metadata which now has the primer info
    writeData(
        wb,
        sheet = "metadata",
        meta_df,
        startRow = 1,
        startCol = 1
    )

    writeData(
        wb,
        sheet = "position_df",
        position_df,
        startRow = 1,
        startCol = 1
    )

    # Save the workbook as an Excel file
    saveWorkbook(
        wb,
        output_file,
        overwrite = TRUE
    )

    print(paste0("Finished! Output file: ", output_file))
}
