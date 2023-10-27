export_biomek_pooling_workbook <- function(assays,
                                           plate_numbers,
                                           position_df_pool,
                                           minipool_overview,
                                           output_dir) {
    # create a blank dataframe to output loop data into
    pooling_df <- data.frame()

    #identify variable for loop below
    minipool_calcs     <- list()
    volume_ranges      <- list()
    minipool_calc_vols <- list()
    minipool_plots     <- list()
    volume_sums        <- c()
    volume_beads       <- c()
    plate_num          <- 0
    out_df             <- data.frame(
        assay = NULL,
        SourcePosition = NULL,
        SourceWell = NULL,
        Volume = NULL,
        DestinationPosition = NULL,
        DestinationWell = NULL
    )

    # for loop that will cycle through assays, plate numbers, and sample types
    for (curr_assay in assays) {
        subs <- position_df_pool |> filter(assay == curr_assay)
        for (curr_plate in sort(plate_numbers)) {
            plate_num <- plate_num + 1
            keys <- c()
            for (sam_type in unique(subs$sample_type)) {
                curr_key   <- paste0(curr_assay, "_", curr_plate, "_", sam_type)
                keys <- append(keys, curr_key)

                # Minipool calculations per plate

                sub_minipool <- minipool_overview |>
                    filter(assay == curr_assay) |>
                    filter(plate_number == curr_plate) |>
                    filter(sample_type == sam_type) |>
                    pull(n_groups)

                minipool_calcs[curr_key] <- subs %>%
                    filter(
                        plate_number == curr_plate &
                        sample_type == sam_type) %>%
                    mutate(miniPool = cut(mean,
                        breaks = sub_minipool),
                miniPool_id = as.integer(cut(mean,
                    breaks = sub_minipool))) %>%
                    data.frame() %>%
                    list()

                # add volume to pool as well as info on tube destination
                # (Plate 1 = P16-A1 well/tube = samples
                # and P16-B1 well/tube = controls,
                # Plate 2 = P16-A2 well/tube = samples
                # and P16-B2 well/tube = controls,)
                volume_ranges[curr_key] <- tibble(
                    miniPool_id = 1:minipool_overview$n_groups[
                        minipool_overview$plate_number == curr_plate &
                        minipool_overview$sample_type == sam_type &
                        minipool_overview$assay == curr_assay
                    ],
                    vol_ul = rev(
                        seq(
                            from = 2,
                            by = 0.5,
                            length.out = minipool_overview$n_groups[
                                minipool_overview$plate_number == curr_plate &
                                minipool_overview$sample_type == sam_type &
                                minipool_overview$assay == curr_assay
                            ]
                        )
                    )
                ) %>%
                    data.frame() %>%
                    list()

                if (sam_type == "sample") {
                    dwell <- paste0("A", curr_plate)
                } else {
                    dwell <- paste0("B", curr_plate)
                }

                biomek_deck_pos <- 12 - plate_num
                sourcepos <- paste0("qPCRPlate-P", biomek_deck_pos)

                minipool_calc_vols[curr_key] <-
                    minipool_calcs[curr_key][[1]] %>%
                    left_join(
                        .,
                        volume_ranges[curr_key][[1]],
                        by = "miniPool_id"
                    ) %>%
                    mutate(
                        DestinationPosition = "P16",
                        DestinationWell = dwell
                    ) %>%
                    data.frame() %>%
                    list()

                # volume_sums  <-
                # append(volume_sums,
                # sum(minipool_calc_vols[curr_key][[1]]$vol_ul))
                # volume_beads <-
                # append(volume_beads,
                # (1.8 * sum(minipool_calc_vols[curr_key][[1]]$vol_ul)))

                # check "minipools" in a plot to confirm group/volume allocation
                minipool_plots[curr_key] <- ggplot(
                        minipool_calc_vols[curr_key][[1]],
                        aes(x = sample_replicate, y = mean)
                    ) +
                    geom_point(aes(colour = as.factor(miniPool_id)), size = 2) +
                    scale_y_continuous(expand = c(0,0)) +
                    xlab("Sample") +
                    ylab("Mean EPF") +
                    theme(axis.text.y = element_text(size = 5)) +
                    # scale_color_brewer("MiniPool ID", palette = "Paired") +
                    coord_flip() %>%
                    list()
            }

            tmp_df <- rbind(
                    minipool_calc_vols[keys[1]][[1]],
                    minipool_calc_vols[keys[2]][[1]]
                ) %>%
                mutate(SourcePosition = sourcepos) %>%
                dplyr::select(
                    assay,
                    SourcePosition,
                    SourceWell = pos,
                    Volume = vol_ul,
                    DestinationPosition,
                    DestinationWell
                ) %>%
                arrange(sample_order(SourceWell))

            out_df <- rbind(out_df, tmp_df)
        }
    }
    write_csv(out_df, paste0(output_dir, "/biomek_pooling_workbook.csv"))

    minipool_all <- minipool_calc_vols %>%
        bind_rows(.id = "MAMAAAAA")
    minipool_all %>%
        mutate(SourcePosition = pos) %>%
        dplyr::select(SourcePosition,
            SourceWell = pos,
            Volume = vol_ul,
            DestinationPosition, DestinationWell) |>
        write_csv(
            paste0(
                output_dir,
                "/biomek_pooling_workbook_mamaaaaa.csv"
            )
        )

    print(paste0("Finished! Output files in: ", output_dir))
}
