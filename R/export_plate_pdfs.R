export_plate_pdfs <- function(df,
                              plate_numbers,
                              assays,
                              output_dir,
                              prefix,
                              plate_width,
                              plate_height) {

    for (curr_plate in plate_numbers) {
        for (curr_assay in assays) {
            pdf(
                paste0(
                    output_dir,
                    "/",
                    prefix,
                    "_",
                    curr_assay,
                    "_",
                    curr_plate,
                    ".pdf"
                )
            )
            plate_plot_epf <- df %>%
                filter(plate_number == curr_plate & assay == curr_assay) %>%
                plate_plot(
                    position = pos,
                    value = EPF,
                    label = round(EPF, digits = 1),
                    plate_size = (plate_width * 2) * (plate_height * 2),
                    title = paste0(unique(.$plate_number), " ", unique(.$assay))
                )
            print(plate_plot_epf)
            dev.off()
        }
    }
}
