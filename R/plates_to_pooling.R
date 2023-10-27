######################################
# Setup
######################################

suppressMessages(library("readxl"))
suppressMessages(library("openxlsx"))
suppressMessages(library("tidyverse"))
suppressMessages(library("dplyr"))
suppressMessages(library("plyr"))
suppressMessages(library("gtools"))
suppressMessages(library("ggplate"))
source("R/sample_order.R")
source("R/read_qPCR_data.R")
source("R/export_plate_pdfs.R")
source("R/export_biomek_pooling_workbook.R")
source("R/print_failed_rep_message.R")
source("R/import_position_df.R")

plates_to_pooling <- function(input_file,
                              qpcr_dir,
                              output_dir,
                              plate_width = 12,
                              plate_height = 8) {
    #### Import LC480 output data ####
    # File name structure is used to determine assay
    # and plate number and so much be YYYYMMDD_Expedition_Assay_Plate#
    # merge all relevant text files from all qPCR plates for a given assay -
    # these will be separate files
    # create a list of file names

    options(dplyr.summarise.inform = FALSE)

    # Make sure directory name doesn't already end in /
    while (endsWith(qpcr_dir, "/")) {
        qpcr_dir <- substr(qpcr_dir, 1, (nchar(qpcr_dir) - 1))
    }

    # Make sure input files exist
    if (!file.exists(qpcr_dir)) {
        stop(
            paste0(
                "Error: ",
                qpcr_dir,
                " doesn't exist."
            )
        )
    }
    if (!file.exists(input_file)) {
        stop(
            paste0(
                "Error: ",
                input_file,
                " doesn't exist."
            )
        )
    }

    # Get all the txt file names in qpcr_dir
    fnames         <- list.files(
        paste0(qpcr_dir, "/"),
        full.names = TRUE,
        pattern = ".txt"
    )

    position_df <- import_position_df(input_file)

    plate_numbers  <- unique(position_df$plate_number)
    assays         <- unique(position_df$assay)

    all_lc480_data <- ldply(fnames, read_qPCR_data, assays, plate_numbers)

    lc480_data_sample <- all_lc480_data %>%
        mutate(plate_number_pos = paste(plate_number, Position, sep = ".")) %>%
        dplyr::select(-Position, -plate_number) %>%
        left_join(., position_df, by = c("plate_number_pos", "assay")) %>%
        filter(., replicate != "pool")

    # TO DO: dont hardcode plate size

    ##### visualise EPF data in a heatmap ####
    prefix <- "raw"
    export_plate_pdfs(
        lc480_data_sample,
        plate_numbers,
        assays,
        output_dir,
        prefix,
        plate_width,
        plate_height
    )

    # identify failed reactions using minimum EPF of 3
    # (for real samples only, criteria not applied to
    # controls as these are all taken through to sequencing )
    rep_failed <- lc480_data_sample %>%
        mutate(
            discard = case_when(
                replicate == "pool" ~ "NA",
                sample_type == "sample" & EPF <= 3 ~ "DISCARD",
                TRUE ~ "KEEP"
            )
        ) %>%
        arrange(sample_order(sample))

    # summary of number of reps to be discarded per sample
    rep_failed_summary <- rep_failed %>%
        filter(replicate != "pool") %>%
        dplyr::group_by(sample, assay) %>%
        dplyr::summarise(count_discard = sum(discard == "DISCARD")) %>%
        ungroup() %>%
        arrange(sample_order(sample))

    # identify the total number of samples per assay
    # with replicates to be discarded
    print_failed_rep_message(rep_failed_summary, assays)

    #identify samples to be completely removed from pool (DISCARD >= 2)
    discarded_samples <- rep_failed_summary %>%
        filter(count_discard >= 2)
    discarded_samples %>%
        arrange(sample)

    # export table for manual removal of failed replicates from 384-well plates

    rep_failed %>%
        dplyr::select(assay, plate_number, pos, sample_replicate, discard) %>%
        filter(discard == "DISCARD") %>%
        arrange(plate_number, sample_order(pos)) %>%
        write_csv(paste0(output_dir, "/reps_to_discard.csv"))

    # visualise plates with failed replicates/samples removed -
    # for sanity check :)

    clean_lc480_data <- rep_failed %>%
        filter(discard == "KEEP")

    ##### visualise clean EPF data in a heatmap ####
    prefix <- "clean"
    export_plate_pdfs(
        clean_lc480_data,
        plate_numbers,
        assays,
        output_dir,
        prefix,
        plate_width,
        plate_height
    )

    # summary EPF samples and controls
    # (at this point is still includes failed replicates -
    # this is ok, as they are filtered later)
    epf_cal <- rep_failed %>%
        dplyr::group_by(assay, sample, sample_type) %>%
        dplyr::summarise(
            median = median(EPF, na.rm = TRUE),
            mean = mean(EPF, na.rm = TRUE),
            sd = sd(EPF, na.rm = TRUE),
            min =  min(EPF),
            max = max(EPF),
            diff_epf = max(EPF) - min(EPF),
            number_valid_reps = sum(EPF >= 3)) %>%
        arrange(sample_order(sample)) %>%
        ungroup()

    # assign calculated mean to respective sample-pool
    # in the original plate layout
    # for samples that have 2 or more useful replicates & all control samples
    epf_cal_mean <- epf_cal %>%
        filter(
            sample_type == "control" |
            (sample_type == "sample" & number_valid_reps >= 2)
        ) %>%
        mutate(sample_replicate = paste0(sample, "-pool")) %>%
        mutate(
            assay_sample_replicate = paste0(
                assay, ".", sample_replicate
            )
        ) %>%
        dplyr::select(assay_sample_replicate, mean)

    position_df_pool <- position_df %>%
        filter(grepl("pool", replicate)) %>%
        mutate(
            assay_sample_replicate = paste0(
                assay, ".", sample_replicate
            )
        ) %>%
        left_join(epf_cal_mean, by = "assay_sample_replicate") %>%
        na.omit(position_df_pool$mean)

    #### Pooling the pooled samples per plate ####
    # fixed volume per plate won't work because of
    # potential big differences in number of samples per plate
    # using grouping of samples based on "minipool" approach (as used in the past)
    # "minipools" created based on 0.5 EPF difference
    # (this is assuming correct florescence values from LC480 where EPF = 15-20 )
    # number of actual samples and controls per plate

    minipool_overview <- position_df_pool %>%
        group_by(assay, plate_number, sample_type) %>%
        dplyr::summarise(
            count_samples = n_distinct(sample),
            min =  min(mean),
            max = max(mean),
            diff_mean = max(mean) - min(mean),
            n_groups = round((max(mean) - min(mean)) / 0.5, digits = 0)
        )

    minipool_overview$n_groups <- ifelse(
        minipool_overview$n_groups < 2,
        2,
        minipool_overview$n_groups
    )

    # Summary of EPF values for actual samples and controls for each plate
    # Where n_groups = the number of minipools each with a range of 0.5 EPF
    # minipool_overview <- position_df_pool %>%
    # group_by(plate_number, sample_type) %>%
    # dplyr::summarise(min =  min(mean), max = max(mean),
    # diff_mean = max(mean) - min(mean),
    # n_groups = round((max(mean) - min(mean))/0.5, digits = 0))
    # minipool_overview

    # Set theme for plots
    theme_set(theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 14, colour = "black"),
            strip.text = element_text(size = 12, colour = "black"),
            legend.position = "bottom",
            legend.text = element_text(size = 11, colour = "black")
        )
    )

    suppressWarnings(
        export_biomek_pooling_workbook(
            assays,
            plate_numbers,
            position_df_pool,
            minipool_overview,
            output_dir
        )
    )

    #library summary with volumes of sample/controls and volume of beads to add for cleanup (1.8 x volume)

    #library_summary %>%

    #cbind(library_summary) %>%
    #    mutate(lib_vol_ul = volume_sums,
    #    bead_vol_ul = volume_beads,
    #    total_vol_ul = lib_vol_ul + bead_vol_ul)

    #confirm total volumes are all <1.5 mL
}
