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

    # We can't forget to add the metadata and position df
    colnames(meta_df)[colnames(meta_df) == 'SAMPLEID'] <- 'sample'
    colnames(meta_df)[colnames(meta_df) == 'SEQUENCINGRUN'] <- 'sequencing_run'
    colnames(meta_df)[colnames(meta_df) == 'ASSAY'] <- 'assay'
    colnames(meta_df)[colnames(meta_df) == 'FW_NO'] <- 'fw_no'
    colnames(meta_df)[colnames(meta_df) == 'RV_NO'] <- 'rv_no'
    colnames(meta_df)[colnames(meta_df) == 'FW_TAG'] <- 'fw_index'
    colnames(meta_df)[colnames(meta_df) == 'RV_TAG'] <- 'rv_index'
    colnames(meta_df)[colnames(meta_df) == 'FW_FULL_SEQ'] <- 'fw_primer'
    colnames(meta_df)[colnames(meta_df) == 'RV_FULL_SEQ'] <- 'rv_primer'
    colnames(meta_df)[colnames(meta_df) == 'PLATE'] <- 'plate'
    colnames(meta_df)[colnames(meta_df) == 'WELL'] <- 'well'
  
    for (i in 1:length(colnames(meta_df))) {
      colnames(meta_df)[i] = tolower(colnames(meta_df)[i])
    }
    meta_df$fastq_1 <- NA
    meta_df$fastq_2 <- NA
    writeData(
        wb,
        sheet = "metadata",
        meta_df,
        startRow = 1,
        startCol = 1
    )
    
    for (assay in assays) {
      addWorksheet(wb, paste0("metadata_", assay))
      
      writeData(
        wb,
        sheet = paste0("metadata_", assay),
        meta_df[meta_df$assay == assay, ],
        startRow = 1,
        startCol = 1
      )
    }

    colnames(position_df) <- c(
        "sample",
        "assay",
        "replicate",
        "pos",
        "sample_replicate",
        "plate_number",
        "assay_plate_number_pos",
        "plate_number_pos",
        "sample_type"
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
