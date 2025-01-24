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
suppressMessages(library("tools"))


##########################################################
# Variables
##########################################################

input_file   <- "test_data/output/output_df.xlsx"
qpcr_dir     <- "test_data/input/qPCR_test_data/"
output_dir   <- "test_data/output/"
plate_width  <- 12
plate_height <- 8
assays       <- c("16S", "MiFish")


##########################################################
# Functions
##########################################################

sample_order <- function(samples) {
  return(order(mixedorder(samples)))
}

read_qPCR_data <- function(file, assays, plate_numbers) {
  extension <- paste0(".", file_ext(basename(file)))
  
  if (extension == ".csv") {
    data <- read.delim(file, sep = ",")
    data <- data.frame(data)
  } else {
    # loads in data
    data <- tryCatch({
      data <- read.delim(file, sep = " ")
      data.frame(data)
    }, error = function(e){
      data <- read.delim(file, sep = "\t")
      data.frame(data)
    }) 
  }

  # get info from sample source, assay and plate id details from file name
  data$fileName <- gsub(
    extension,
    "",
    basename(file)
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

  
  if (! curr_plate %in% plate_numbers[curr_assay][[1]]) {
    stop(
      paste0(
        "Error: plate - ",
        curr_plate,
        " taken from file ",
        file,
        " not found in ",
        plate_numbers[curr_assay],
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

  for (curr_assay in assays) {
    for (curr_plate in plate_numbers[curr_assay][[1]]) {
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
          value = epf,
          label = round(epf, digits = 1),
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

  # identify variables for loop below
  minipool_calcs     <- list()
  volume_ranges      <- list()
  minipool_calc_vols <- list()
  minipool_plots     <- list()
  volume_sums        <- c()
  volume_beads       <- c()
  plate_num          <- 0
  book               <- 0
  sam_type_num       <- 0
  biomek_deck_pos    <- 1
  assay_num          <- 0
  out_df             <- data.frame(
    assay = NULL,
    plate_number = NULL,
    SourcePosition = NULL,
    SourceWell = NULL,
    Volume = NULL,
    DestinationPosition = NULL,
    DestinationWell = NULL
  )
  out_dfs            <- list()

  # Loop that will cycle through assays, plate numbers, and sample types
  for (curr_assay in assays) {
    subs      <- position_df_pool |> filter(assay == curr_assay)
    assay_num <- assay_num + 1
    for (curr_plate in plate_numbers[curr_assay][[1]]) {
      unique_sam_types <- unique(subs$sample_type)
      plate_num        <- plate_num + 1
      keys             <- c()

      for (sam_type in unique_sam_types) {
        curr_key <- paste0(curr_assay, "_", curr_plate, "_", sam_type)
        keys     <- append(keys, curr_key)

        # Minipool calculations per plate per assay

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

        volume_ranges[curr_key] <- tibble(
          miniPool_id = 1:6,
          vol_ul = rev(
            seq(
              from = 4,
              by = 1,
              length.out = 6
            )
          )
        ) %>%
          data.frame() %>%
          list()
        #  p_num <- substr(curr_plate, 6, nchar(curr_plate))


        sam_type_num <- sam_type_num + 1
        if (sam_type_num > length(unique_sam_types)) {
          sam_type_num    <- 1
          biomek_deck_pos <- biomek_deck_pos + 1
          if (biomek_deck_pos > 4) {
            biomek_deck_pos    <- 1
            book               <- book + 1
            out_dfs[[book]]    <- out_df

            out_df             <- data.frame(
              assay = NULL,
              plate_number = NULL,
              SourcePosition = NULL,
              SourceWell = NULL,
              Volume = NULL,
              DestinationPosition = NULL,
              DestinationWell = NULL
            )
          }
        }
        sourcepos <- paste0("qPCR", biomek_deck_pos)

        if (sam_type == "sample") {
          dwell <- paste0("A", biomek_deck_pos)
        } else {
          dwell <- paste0("B", biomek_deck_pos)
        }

        # This try catch exists for when a plate is missing a sampletype
        minipool_calc_vols[curr_key] <- tryCatch({
          minipool_calcs[curr_key][[1]] %>%
            left_join(
              .,
              volume_ranges[curr_key][[1]],
              by = "miniPool_id"
            ) %>%
            mutate(
              DestinationPosition = "PoolTubes",
              DestinationWell = dwell
            ) %>%
            data.frame() %>%
            list()
        }, error = function(e){
          minipool_calcs[curr_key][[1]] %>%
            left_join(
              .,
              volume_ranges[curr_key][[1]],
              by = "miniPool_id"
            ) %>%
            mutate(
              DestinationPosition = NULL,
              DestinationWell = NULL
            ) %>%
            data.frame() %>%
            list()
        })

        minipool_calc_vols[curr_key][[1]] <- minipool_calc_vols[curr_key][[1]] %>%
          mutate(DestinationWell = ifelse(grepl("^ITC_", SAMPLE), paste0("C", assay_num), DestinationWell))

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
        #  check "minipools" in a plot to confirm group/volume allocation
        minipool_plots[[curr_key]] <- ggplot(
          minipool_calc_vols[curr_key][[1]],
          aes(x = sample_replicate, y = mean)
        ) +
          geom_point(aes(colour = as.factor(miniPool_id)), size = 2) +
          scale_y_continuous(expand = c(0,0)) +
          xlab("Sample") +
          ylab("Mean EPF") +
          theme(axis.text.y = element_text(size = 5)) +
          coord_flip()

        if (sam_type_num == length(unique_sam_types)) {
          if (length(unique_sam_types) == 1) {
            tmp_df <- minipool_calc_vols[keys[1]][[1]] %>%
            mutate(SourcePosition = sourcepos) %>%
            dplyr::select(
              assay,
              plate_number,
              SourcePosition,
              SourceWell = pos,
              Volume = vol_ul,
              DestinationPosition,
              DestinationWell
            ) %>%
            arrange(sample_order(SourceWell))
          } else if (length(unique_sam_types) == 2) {
            tmp_df <- rbind(
              minipool_calc_vols[keys[1]][[1]],
              minipool_calc_vols[keys[2]][[1]]
            ) %>%
            mutate(SourcePosition = sourcepos) %>%
            dplyr::select(
              assay,
              plate_number,
              SourcePosition,
              SourceWell = pos,
              Volume = vol_ul,
              DestinationPosition,
              DestinationWell
            ) %>%
            arrange(sample_order(SourceWell))
          }

          out_df <- rbind(out_df, tmp_df)
        }

        pdf(paste0(output_dir, curr_key, "_minipool_plot.pdf"))
        print(minipool_plots[curr_key])
        dev.off()
      }
    }
  }
  book              <- book + 1
  out_dfs[[book]]   <- out_df

  for (i in 1:book) {
    write_csv(out_dfs[[i]], paste0(output_dir, "/biomek_pooling_workbook_", i,".csv")) 
  }

  print(paste0("Finished! Output files in: ", output_dir))

  return(minipool_calc_vols)
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
  cat("--------------------------------------------------------------------------------")
  cat("\n\n")
}

change_lc480_colnames <- function(lc480) {
  colnames(lc480) <- toupper(colnames(lc480))
  colnames(lc480)[which(names(lc480) == "EXPERIMENT_NAME")] <- "experiment_name"
  colnames(lc480)[which(names(lc480) == "POSITION")] <- "position"
  colnames(lc480)[which(names(lc480) == "WELL_POSITION")] <- "position"
  colnames(lc480)[which(names(lc480) == "SAMPLE")] <- "sample"
  colnames(lc480)[which(names(lc480) == "EPF")] <- "epf"
  colnames(lc480)[which(names(lc480) == "DRN")] <- "epf"
  colnames(lc480)[which(names(lc480) == "CP")] <- "cp"
  colnames(lc480)[which(names(lc480) == "CQ")] <- "cp"
  colnames(lc480)[which(names(lc480) == "TM1")] <- "tm1"
  colnames(lc480)[which(names(lc480) == "TM2")] <- "tm2"
  colnames(lc480)[which(names(lc480) == "FILENAME")] <- "filename"
  colnames(lc480)[which(names(lc480) == "ASSAY")] <- "assay"
  colnames(lc480)[which(names(lc480) == "PLATE_NUMBER")] <- "plate_number"
  return(lc480)
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
  colnames(position_df) <- toupper(colnames(position_df))

  # Validata columns
  if (!("ASSAY" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain an 'assay' column")
  } else {
    colnames(position_df)[which(names(position_df) == "ASSAY")] <- "assay"
  }
  if (!("PLATE_NUMBER" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain a 'plate_number' column")
  } else {
    colnames(position_df)[which(names(position_df) == "PLATE_NUMBER")] <- "plate_number"
  }
  if (!("PLATE_NUMBER_POS" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain a 'plate_number_pos' column")
  } else {
    colnames(position_df)[which(names(position_df) == "PLATE_NUMBER_POS")] <- "plate_number_pos"
  }
  if (!("REPLICATE" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain a 'replicate' column")
  } else {
    colnames(position_df)[which(names(position_df) == "REPLICATE")] <- "replicate"
  }
  if (!("SAMPLE_REPLICATE" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain a 'sample_replicate' column")
  } else {
    colnames(position_df)[which(names(position_df) == "SAMPLE_REPLICATE")] <- "sample_replicate"
  }
  if (!("SAMPLE_TYPE" %in% colnames(position_df))) {
    stop("Sheet 'position_df' must contain a 'sample_type' column")
  } else {
    colnames(position_df)[which(names(position_df) == "SAMPLE_TYPE")] <- "sample_type"
  }
  colnames(position_df)[which(names(position_df) == "SAMPLE_ID")] <- "sample_id"
  colnames(position_df)[which(names(position_df) == "POS")] <- "pos"
  colnames(position_df)[which(names(position_df) == "ASSAY_PLATE_NUMBER_POS")] <- "assay_plate_number_pos"

  return(position_df)
}

import_samplesheet_df <- function(excel_file, assay) {
  excel_sheets   <- excel_sheets(excel_file)
  metadata_sheet <- grep(
    glob2rx(paste0("metadata_", assay)),
    excel_sheets,
    ignore.case = TRUE,
    value = TRUE
  )
  
  # Validate number of sheets called 'position_df'
  if (length(metadata_sheet) == 0) {
    stop(paste0("No sheets found named 'metadata_'", assay, "'"))
  }
  if (length(metadata_sheet) > 1) {
    stop(paste0("Multiple sheets found named 'metadata_'", assay, "'"))
  }
  
  # Import df
  metadata_df   <- read_excel(excel_file, sheet = paste0("metadata_", assay)) %>%
    as.data.frame()
  
  return(metadata_df)
}


##########################################################
# Import data
##########################################################

if (TRUE) {  
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
    pattern = ".txt|.csv|.tsv"
  )

  fnames         <- str_replace(fnames, "//", "/")
  position_df    <- import_position_df(input_file)
  
  plate_numbers <- list()
  for (assay in assays) {
    plate_numbers[assay] = list(unique(position_df[position_df$assay == assay, "plate_number"]))
  }

  all_lc480_data <- ldply(fnames, read_qPCR_data, assays, plate_numbers)
  all_lc480_data <- change_lc480_colnames(all_lc480_data)
  
  all_lc480_data$sample <- tryCatch({
    all_lc480_data$sample <- str_replace(all_lc480_data$sample, "Sample", "")
    all_lc480_data$sample <- trimws(all_lc480_data$sample)
    all_lc480_data$sample
  }, error = function(e) {
    all_lc480_data$sample <- seq_len(nrow(all_lc480_data))
    all_lc480_data$sample
  })
  
  # Make sure qPCR data isn't missing any plates
  for (assay in assays) {
    missing_plates <- c()
    for (plate in plate_numbers[assay][[1]]) {
      if (! plate %in% unique(all_lc480_data[all_lc480_data$assay == assay, "plate_number"])) {
        missing_plates <- c(missing_plates, plate)
      }
    }
    if (length(missing_plates) > 0) {
      print("qPCR data directory is missing plates: ")
      print(missing_plates)
      stop()
    }
  }

  lc480_data_sample <- all_lc480_data %>%
    mutate(plate_number_pos = paste(plate_number, position, sep = ".")) %>%
    dplyr::select(-position, -plate_number) %>%
    left_join(., position_df, by = c("plate_number_pos", "assay")) %>%
    filter(., replicate != "pool")

  lc480_data_sample$sample <- lc480_data_sample$SAMPLE
}


##########################################################
# Visualise raw data
##########################################################

# visualise EPF data in a heatmap
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

# The plots are outputed into the output_dir specified as .pdf but let's look at one of those plots in R studio. 
plate_plot_epf <- lc480_data_sample %>%
  filter(plate_number == "Plate2" & assay == "16S") %>%
  plate_plot(
    position = pos,
    value = epf,
    label = round(epf, digits = 1),
    plate_size = (plate_width * 2) * (plate_height * 2),
    title = paste0(unique(.$plate_number), " ", unique(.$assay))
  )
print(plate_plot_epf)


##########################################################
# Flag samples for discarding or double checking
##########################################################

# identify failed reactions using minimum EPF of 2 and maximum cp of 40
# (for real samples only, criteria not applied to
# controls as these are all taken through to sequencing )
rep_failed <- lc480_data_sample %>%
  mutate(
    discard = case_when(
      replicate == "pool" ~ "NA",
      (sample_type == "sample" & epf < 2) | (sample_type == "sample" & cp > 40) ~ "DISCARD",
      TRUE ~ "KEEP"
    )
  ) %>%
  arrange(sample_order(sample))

# Flag samples with epf between 2 and 5
rep_failed <- rep_failed %>%
  mutate(
    epf_2_to_5 = case_when(
      (sample_type == "sample" & epf >= 2) & (sample_type == "sample" & epf <= 5) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  arrange(sample_order(sample))

# Flag samples with tm1 outside of standard deviation
rep_failed$tm1_sd         <- NA
rep_failed$tm1_mean       <- NA
rep_failed$tm1_outside_sd <- NA
for (curr_assay in assays) {
  plates <- unique(rep_failed[rep_failed$assay == curr_assay, "plate_number"])
  for (curr_plate in plates) {
    tm1s       <- rep_failed[rep_failed$assay == curr_assay & rep_failed$plate_number == curr_plate & rep_failed$sample_type == "sample", "tm1"]
    tm1_sd     <- sd(na.omit(tm1s))
    tm1_mean   <- mean(na.omit(tm1s))
    tm1_min    <- tm1_mean - tm1_sd
    tm1_max    <- tm1_mean + tm1_sd
    rep_failed[rep_failed$assay == curr_assay & rep_failed$plate_number == curr_plate, "tm1_sd"]   <- sd(na.omit(tm1s))
    rep_failed[rep_failed$assay == curr_assay & rep_failed$plate_number == curr_plate, "tm1_mean"] <- mean(na.omit(tm1s))

    rep_failed[rep_failed$assay == curr_assay & rep_failed$plate_number == curr_plate, ] <- rep_failed[rep_failed$assay == curr_assay & rep_failed$plate_number == curr_plate, ] %>%
      mutate(
        tm1_outside_sd = case_when(
          (tm1 < tm1_min | tm1 > tm1_max) ~ TRUE,
          TRUE ~ FALSE
        )
      ) %>%
      arrange(sample_order(sample))
  }
}


##########################################################
# View the worse reps
##########################################################

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

# print the failed reps (reps where the epf < 2 and the Cp >40)
discard_df <- print(rep_failed[rep_failed$discard == "DISCARD", ])

# print the reps with epf between 2 to 5  
print(rep_failed[rep_failed$epf_2_to_5 == TRUE, ])

# print reps with tm1 outside of standard deviation
print(rep_failed[rep_failed$tm1_outside_sd == TRUE, ])

# Print all the bad reps
print(rep_failed[(rep_failed$tm1_outside_sd == TRUE | rep_failed$epf_2_to_5 == TRUE | rep_failed$discard == "DISCARD") & rep_failed$sample_type == "sample", ])


##########################################################
# Output rxns_to_check.csv
##########################################################

write.csv(row.names = FALSE, rep_failed[(rep_failed$tm1_outside_sd == TRUE | rep_failed$epf_2_to_5 == TRUE | rep_failed$discard == "DISCARD") & rep_failed$sample_type == "sample", ], paste0(output_dir, "rxns_to_check.csv"))


######################################################
# Import rxns_to_check.csv and merge discard column
######################################################

checked_runs <- read.csv(paste0(output_dir, "rxns_to_check.csv"))

for (row in 1:nrow(checked_runs)) {
  id <- checked_runs[row, "assay_plate_number_pos"]
  disc <- toupper(checked_runs[row, "discard"])

  if (disc != "KEEP" & disc != "DISCARD") {
    while (TRUE) {
      print(paste0("discard value of ", disc, " is not valid for ", id))
      answer <- readline(prompt="Enter 1 to keep, or 2 to discard: ")
      if (answer == 1) {
        disc <- "KEEP"
        break
      } else if (answer == 2) {
        disc <- "DISCARD"
        break
      } else {
        print(paste0(answer, " is not a valid option"))
      }
    }
  }

  rep_failed$discard[rep_failed$assay_plate_number_pos == id] <- disc
}

write.csv(row.names = FALSE, rep_failed, paste0(output_dir, "rxns_checked.csv"))

# summary of number of reps to be discarded per sample 
rep_failed_summary <- rep_failed %>%   
  filter(replicate != "pool") %>%   
  dplyr::group_by(sample, assay) %>%   
  dplyr::summarise(count_discard = sum(discard == "DISCARD")) %>%   
  ungroup() %>%   
  arrange(sample_order(sample)) 
# identify the total number of samples per assay 
# with replicates to be discarded print_failed_rep_message(rep_failed_summary, assays) 
# print the failed reps (reps where the epf < 2 and the Cp >40) 
discard_df <- print(rep_failed[rep_failed$discard == "DISCARD", ]) 
# print the reps with epf between 2 to 5   
print(rep_failed[rep_failed$epf_2_to_5 == TRUE, ]) 
# print reps with tm1 outside of standard deviation 
print(rep_failed[rep_failed$tm1_outside_sd == TRUE, ]) 
# Print all the bad reps 
print(rep_failed[(rep_failed$tm1_outside_sd == TRUE | rep_failed$epf_2_to_5 == TRUE | rep_failed$discard == "DISCARD") & rep_failed$sample_type == "sample", ]) 


##########################################################
# Output reps_to_discard.csv
##########################################################

#identify samples to be completely removed from pool (DISCARD >= 2)
discarded_samples <- rep_failed_summary %>%
  filter(count_discard >= 2) %>%
  arrange(sample)

# export table for manual removal of failed replicates from 384-well plates
## AB - can this please only include samples where 1 replicate needs to be discarded.
# In cases where >1 rep is discarded
reps_to_discard <- rep_failed %>%
  dplyr::select(assay, plate_number, pos, sample_replicate, discard) %>%
  filter(discard == "DISCARD") %>%
  arrange(plate_number, sample_order(pos))

reps_to_discard$sample <- substr(reps_to_discard$sample_replicate, 1, nchar(reps_to_discard$sample_replicate)-2)

# Keep only rows where the sample is not duplicated
duplicated_rows <- duplicated(reps_to_discard[c("assay", "sample")]) | duplicated(reps_to_discard[c("assay", "sample")], fromLast = TRUE)
reps_to_discard <- subset(reps_to_discard, !duplicated_rows)

write_csv(reps_to_discard, paste0(output_dir, "/reps_to_discard.csv"))


##########################################################
# Visualise clean data
##########################################################

clean_lc480_data <- rep_failed %>%
  dplyr::group_by(assay, sample) %>%
  filter(discard=="KEEP")


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

# Let's look at one of those plots in R studio
clean_plate_plot_epf <- clean_lc480_data %>%
  filter(plate_number == "Plate5" & assay == "MiFish") %>%
  plate_plot(
    position = pos,
    value = epf,
    label = round(epf, digits = 1),
    plate_size = (plate_width * 2) * (plate_height * 2),
    title = paste0(unique(.$plate_number), " ", unique(.$assay))
  )
print(clean_plate_plot_epf)


##########################################################
# EPF calcs
##########################################################

# summary EPF samples and controls
# (at this point is still includes failed replicates -
# this is ok, as they are filtered later)
epf_cal <- rep_failed %>%
  dplyr::group_by(assay, sample, sample_type) %>%
  dplyr::summarise(
    median = median(epf, na.rm = TRUE),
    mean = mean(epf, na.rm = TRUE),
    sd = sd(epf, na.rm = TRUE),
    min =  min(epf),
    max = max(epf),
    diff_epf = max(epf) - min(epf),
    number_valid_reps = (sum(discard == "KEEP") >= 2)) %>%
  arrange(sample_order(sample)) %>%
  ungroup()

# assign calculated mean to respective sample-pool
# in the original plate layout
# for samples that have 2 or more useful replicates & all control samples
epf_cal_mean <- epf_cal %>%
  filter(
    sample_type == "control" |
      (sample_type == "sample" & number_valid_reps == TRUE)
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
  left_join(epf_cal_mean, by = "assay_sample_replicate") #%>%
  #na.omit(position_df_pool$mean)
position_df_pool <- position_df_pool[(position_df_pool$sample_type == "control" | !is.na(position_df_pool$mean)), ]

##########################################################
# Minipool Generation Calcs & Run output function
##########################################################

# These calculations use grouping of samples based on "minipool" approach 
# For each plate 6 minipools are created for the samples

minipool_overview <-   position_df_pool %>%
  group_by(assay, plate_number, sample_type) %>%
  dplyr::summarise(
    count_samples = n_distinct(SAMPLE),
    min =  min(mean),
    max = max(mean),
    diff_mean = max(mean) - min(mean),
    n_groups = 6
  )
position_df_pool_filt <- position_df_pool[!is.na(position_df_pool$mean),]
minipool_overview_filt <- minipool_overview[!is.na(minipool_overview$diff_mean),] 
minipool_calc_vols <- export_biomek_pooling_workbook(assays, plate_numbers, position_df_pool_filt, minipool_overview_filt, output_dir)

# Minipool summaries including final volume in pooled tube & the number of samples/controls in each tube
minipool_vols_df      <- do.call(rbind.data.frame, minipool_calc_vols)
minipool_vols_summary <- minipool_vols_df %>%
  group_by(assay, plate_number, sample_type, DestinationWell) %>%
  dplyr::summarise(
    total_vol_ul = sum(vol_ul),
    count_samples = n_distinct(SAMPLE)) %>%
  ungroup() %>%
  mutate(assay.plate_number.sample_type = paste(assay, plate_number, sample_type, sep = ".") )

write_csv(minipool_vols_summary, paste0(output_dir, "/minipool_summary.csv"))

#confirm final tube volumes do not exceed 1.5 mL (or 1500 uL)


# Create samplesheets with info on discarded samples
for (assay in assays) {
  meta_df           <- import_samplesheet_df(input_file, assay)
  meta_df$discarded <- FALSE
  curr_disc_sams    <- discarded_samples[discarded_samples$assay == assay, ]
  
  meta_df$discarded <- ifelse(meta_df$sample %in% curr_disc_sams$sample, TRUE, meta_df$discarded)
  
  write_csv(meta_df, paste0(output_dir, "/samplesheet_", assay, ".csv"))
}
