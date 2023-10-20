######################################
# Setup
######################################

library("readxl")
library("openxlsx")
library("tidyverse")
library("dplyr")
library("plyr")
library("gtools")
library("ggplate")
#source("R/sample_order.R")
#source("R/read_qPCR_data.R")
sample_order <- function(samples) {
    return(order(mixedorder(samples)))
}

# Loop to manipulate all files simultaneously
read_qPCR_data <- function(file) {
    # loads in data
    data <- read.delim(file)

    # get info from sample source, assay and plate id details from file name
    data$fileName <- gsub(
        ".txt",
        "",
        str_split(file, "\\/\\/", simplify = TRUE)[, 2]
    )
    description       <- strsplit(data$fileName, "_")
    data$date         <- sapply(description, "[", 1)
    data$expedition   <- sapply(description, "[", 2)
    data$assay        <- sapply(description, "[", 3)
    data$plate_number <- sapply(description, "[", 4)

    return(data)
}


#### Import LC480 output data ####
# File name structure is used to determine assay
# and plate number and so much be YYYYMMDD_Expedition_Assay_Plate#
# merge all relevant text files from all qPCR plates for a given assay -
# these will be separate files
# create a list of file names
exp_name <- "testing_code"
fnames   <- list.files("qPCR_data/", full.names = TRUE, pattern = ".txt")

all_lc480_data <- ldply(fnames, read_data)
position_df    <- read_excel("results.xlsx", sheet = "position_df") %>%
    as.data.frame()

lc480_data_sample <- all_lc480_data %>%
    mutate(plate_number_pos = paste(plate_number, Position, sep = ".")) %>%
    dplyr::select(-Position, -plate_number) %>%
    left_join(., position_df, by = c("plate_number_pos", "assay")) %>%
    filter(., replicate != "pool")

##### visualise EPF data in a heatmap ####
# You will need to change the plate_number and assay for your relevant dataset
# Plate 1 16S

plate1_16S_epf <- lc480_data_sample %>%
    filter(plate_number == "Plate1" & assay == "16S") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay))
    )

#Plate 2 16S
plate2_16S_epf <- lc480_data_sample %>%
    filter(plate_number == "Plate2" & assay == "16S") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay))
    )

# Plate 1 MiFish

plate1_MiFish_epf <- lc480_data_sample %>%
    filter(plate_number == "Plate1" & assay == "MiFish") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay))
    )

#Plate 2 MiFish
plate2_MiFish_epf <- lc480_data_sample %>%
    filter(plate_number == "Plate2" & assay == "MiFish") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay))
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

# identify the total number of samples per assay with replicates to be discarded 
#16S
cat(
    "There are ",
    (rep_failed_summary %>%
        filter(assay == "16S", count_discard == 1) %>%
        n_distinct()
    ),
    " samples in 16S assay with 1 replicate to be discarded\n",
    "There are ",
    (rep_failed_summary %>%
        filter(assay == "16S", count_discard >= 2) %>%
        n_distinct()),
    " samples in 16S assay with 2 or more failed replicates\n"
)

#MiFish
cat("There are ",
    (rep_failed_summary %>%
        filter(assay == "MiFish", count_discard == 1) %>%
        n_distinct()),
    " samples in MiFish assay with 1 replicate to be discarded\n",
    "There are ",
    (rep_failed_summary %>%
        filter(assay == "MiFish", count_discard >= 2) %>%
        n_distinct()),
    " samples in MiFish assay with 2 or more failed replicates\n"
)


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
    write_csv("reps_to_discard.csv")

# visualise plates with failed replicates/samples removed - for sanity check :)

clean_lc480_data <- rep_failed %>%
    filter(discard == "KEEP")

#Plate 1 16S Clean
plate1_16S_clean <- clean_lc480_data %>%
    filter(plate_number == "Plate1", assay == "16S") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay), " Clean"),
    )



#Plate 2 16S clean
plate2_16S_clean <- clean_lc480_data %>%
    filter(plate_number == "Plate2", assay == "16S") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay), " Clean"),
    )


#Plate 1 MiFish clean
plate1_MiFish_clean <- clean_lc480_data %>%
    filter(plate_number == "Plate1", assay == "MiFish") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay), " Clean"),
    )

#Plate 2 MiFish clean
plate2_MiFish_clean <- clean_lc480_data %>%
    filter(plate_number == "Plate2", assay == "MiFish") %>%
    plate_plot(
        position = pos,
        value = EPF,
        label = round(EPF, digits = 1),
        plate_size = 384,
        title = paste0(unique(.$plate_number), " ", unique(.$assay), " Clean"),
    )

# Instead of the old way, lets just produce pdfs of all the plots
#clean_lc480_data %>%
    #filter(plate_number == "Plate2", assay == "MiFish") %>%
    #plate_pdf(plate_width, plate_height, paste0("Plate2", "_plot_", "MiFish", ".pdf"))





# summary EPF samples and controls (at this point is still includes failed replicates - this is ok, as they are filtered later)
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

# assign calculated mean to respective sample-pool in the original plate layout 
# for samples that have 2 or more useful replicates & all control samples
epf_cal_mean <- epf_cal %>%
    filter(
        sample_type == "control" |
        (sample_type == "sample" & number_valid_reps >= 2)
    ) %>%
    mutate(sample_replicate = paste0(sample, "-pool")) %>%
    mutate(assay_sample_replicate = paste0(assay, ".", sample_replicate)) %>%
    dplyr::select(assay_sample_replicate, mean)

position_df_pool <- position_df %>%
    filter(grepl("pool", replicate)) %>%
    mutate(assay_sample_replicate = paste0(assay, ".", sample_replicate)) %>%
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
        n_groups = round((max(mean) - min(mean))/0.5, digits = 0)
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
    assay=NULL,
    SourcePosition=NULL,
    SourceWell=NULL,
    Volume=NULL,
    DestinationPosition=NULL,
    DestinationWell=NULL
)

# for loop that will cycle through assays, plate numbers, and sample types
for (curr_assay in unique(position_df_pool$assay)) {
    subs <- position_df_pool |> filter(assay == curr_assay)
    for (plate in unique(subs$plate_number)) {
        plate_num  <- plate_num + 1
        plate_keys <- c()
        for (sam_type in unique(subs$sample_type)) {
            curr_key   <- paste0(curr_assay, "_", plate, "_", sam_type)
            plate_keys <- append(plate_keys, curr_key)

            # Minipool calculations per plate

            sub_minipool <- minipool_overview |>
                filter(assay == curr_assay) |>
                filter(plate_number == plate) |>
                filter(sample_type == sam_type) |>
                pull(n_groups)

            minipool_calcs[curr_key] <- subs %>%
                filter(plate_number == plate & sample_type == sam_type) %>%
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
                    minipool_overview$plate_number == plate &
                    minipool_overview$sample_type == sam_type &
                    minipool_overview$assay == curr_assay
                ],
                vol_ul = rev(
                    seq(
                        from = 2,
                        by = 0.5,
                        length.out = minipool_overview$n_groups[
                            minipool_overview$plate_number == plate &
                            minipool_overview$sample_type == sam_type &
                            minipool_overview$assay == curr_assay
                        ]
                    )
                )
            ) %>%
                data.frame() %>%
                list()

            if (sam_type == "sample") {
                dwell <- paste0("A", plate_num)
            } else {
                dwell <- paste0("B", plate_num)
            }

            biomek_deck_pos <- 12 - plate_num
            sourcepos <- paste0("qPCRPlate-P", biomek_deck_pos)

            minipool_calc_vols[curr_key] <- minipool_calcs[curr_key][[1]] %>%
                left_join(
                    .,
                    volume_ranges[curr_key][[1]],
                    by = "miniPool_id"
                ) %>%
                mutate(DestinationPosition = "P16", DestinationWell = dwell) %>%
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
                minipool_calc_vols[plate_keys[1]][[1]],
                minipool_calc_vols[plate_keys[2]][[1]]
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
            arrange(sample_order(SourceWell)) #%>%
            #write_csv(paste0("biomek_pooling_workbook_", curr_assay, plate,".csv"))

        out_df <- rbind(out_df, tmp_df)
    }

}
write_csv(out_df, "biomek_pooling_workbook.csv")

minipool_all <- minipool_calc_vols %>%
    bind_rows(.id = 'MAMAAAAA') 
minipool_all %>%
    mutate(SourcePosition = pos) %>%
    dplyr::select(SourcePosition,
        SourceWell = pos, 
        Volume = vol_ul, 
        DestinationPosition, DestinationWell) |> 
    write_csv(paste0("biomek_pooling_workbook_", exp_name,".csv"))
    #View()

#minipool_plots["16S_Plate1_sample"][[1]]
#minipool_plots["16S_Plate1_control"][[1]]
#minipool_plots["COI_Plate1_sample"][[1]]
#minipool_plots["COI_Plate1_control"][[1]]


#library summary with volumes of sample/controls and volume of beads to add for cleanup (1.8 x volume)

#library_summary %>%

cbind(library_summary) %>%
    mutate(lib_vol_ul = volume_sums,
    bead_vol_ul = volume_beads,
    total_vol_ul = lib_vol_ul + bead_vol_ul)

#confirm total volumes are all <1.5 mL
