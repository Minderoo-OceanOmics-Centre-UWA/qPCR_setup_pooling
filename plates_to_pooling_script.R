##########################################################
# Packages
##########################################################

suppressMessages(library("readxl"))
suppressMessages(library("openxlsx"))
suppressMessages(library("tidyverse"))
suppressMessages(library("dplyr"))
suppressMessages(library("plyr"))
suppressMessages(library("gtools"))
suppressMessages(library("ggplate"))


##########################################################
# Variables
##########################################################

input_file   <- "test_data/results.xlsx"
qpcr_dir     <- "test_data/qPCR_test_data"
output_dir   <- "test_data/pooling_dir"
plate_width  <- 12
plate_height <- 8


##########################################################
# Functions
##########################################################

sample_order <- function(samples) {
  return(order(mixedorder(samples)))
}

read_qPCR_data <- function(file, assays, plate_numbers) {
  # loads in data
  data <- read.delim(file)
  
  # get info from sample source, assay and plate id details from file name
  data$fileName <- gsub(
    ".txt",
    "",
    str_split(file, "\\/\\/", simplify = TRUE)[, 2]
  )
  description       <- strsplit(data$fileName, "_")
  desc_count        <- length(description[[1]])
  data$assay        <- sapply(description, "[", (desc_count - 1))
  data$plate_number <- sapply(description, "[", desc_count)
  
  curr_assay <- unique(data$assay)
  curr_plate <- unique(data$plate_number)
  if (! curr_assay %in% assays) {
    stop(
      paste0(
        "Error: assay - ",
        curr_assay,
        " taken from file ",
        file,
        " not found in ",
        assays,
        " Make sure filename ends in _$assay_$plate.txt"
      )
    )
  }
  
  if (! curr_plate %in% plate_numbers) {
    stop(
      paste0(
        "Error: plate - ",
        curr_plate,
        " taken from file ",
        file,
        " not found in ",
        plate_numbers,
        " Make sure filename ends in _$assay_$plate.txt"
      )
    )
  }
  
  return(data)
}

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

print_failed_rep_message <- function(rep_failed_summary, assays) {
  cat("\n")
  cat("--------------------------------------------------------------------------------")
  cat("\nFailed Replicate Message: ")
  for (curr_assay in assays) {
    count_eq_one <- rep_failed_summary %>%
      filter(assay == curr_assay, count_discard == 1) %>%
      n_distinct()
    count_ge_two <- rep_failed_summary %>%
      filter(assay == curr_assay, count_discard >= 2) %>%
      n_distinct()
    cat(
      (ifelse(count_eq_one == 1, "\nThere is", "\nThere are")),
      count_eq_one,
      (ifelse(count_eq_one == 1, "sample in", "samples in")),
      curr_assay,
      "assay with 1 replicate to be discarded",
      (ifelse(count_ge_two == 1, "\nThere is", "\nThere are")),
      count_ge_two,
      (ifelse(count_ge_two == 1, "sample in", "samples in")),
      curr_assay,
      "assay with 2 or more failed replicates\n"
    )
  }
  #cat("\n")
  cat("--------------------------------------------------------------------------------")
  cat("\n\n")
}

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


##########################################################
# Import data
##########################################################
  
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
fnames <- list.files(
  paste0(qpcr_dir, "/"),
  full.names = TRUE,
  pattern = ".txt"
)
  
position_df <- import_position_df(input_file)


##########################################################
# Main
##########################################################
  
plate_numbers  <- unique(position_df$plate_number)
assays         <- unique(position_df$assay)
  
all_lc480_data <- ldply(fnames, read_qPCR_data, assays, plate_numbers)

# Make sure qPCR data isn't missing any plates
missing_plates <- c()
for (plate in plate_numbers) {
  if (! plate %in% unique(all_lc480_data$plate_number)) {
    missing_plates <- c(missing_plates, plate)
  }
}
if (length(missing_plates) > 0) {
  print("qPCR data directory is missing plates: ")
  print(missing_plates)
  stop()
}
  
lc480_data_sample <- all_lc480_data %>%
  mutate(plate_number_pos = paste(plate_number, Position, sep = ".")) %>%
  dplyr::select(-Position, -plate_number) %>%
  left_join(., position_df, by = c("plate_number_pos", "assay")) %>%
  filter(., replicate != "pool")
  
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

