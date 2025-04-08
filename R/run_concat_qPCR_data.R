##########################################################
# Libraries
##########################################################

library(readr)
library(stringr)


##########################################################
# Variables
##########################################################

input_dir   <- "test_data/input/cp_epf_tm_files/"
output_dir  <- "test_data/output/"
assays      <- c("16S", "MiFish")
start_plate <- 1


##########################################################
# Main
##########################################################

for (assay in assays) {
  cp_files    <- Sys.glob(paste0(input_dir, "*", assay, "*Cp.txt"))
  plate_count <- length(cp_files)
  
  for (plate_num in start_plate:(start_plate + plate_count - 1)) {
    curr_plate <- paste0("Plate", plate_num)
    print(paste0(assay, ": ", curr_plate))
    
    cp_file  <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*Cp.txt"))
    epf_file <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*EPF.txt"))
    tm_file  <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*Tm.txt"))
    
    # Import cp and keep only Position, and Cp
    print(paste0("Importing: ", cp_file))
    cp          <- read_delim(cp_file, delim = "\t", skip = 1)
    if (dim(cp)[1] == 0) {
      stop(paste0("cp data appears to be empty. Check that this path is correct and the file isn't empty ", cp_file))
    }
    cp$Position <- cp$Pos
    cp          <- cp[, c("Position", "Cp")]
    
    # Import epf, remove weird columns, transpose, keep final EPF plus Position and Sample
    print(paste0("Importing: ", epf_file))
    epf                          <- read_delim(epf_file, delim = "\t")
    if (dim(epf)[1] == 0) {
      stop(paste0("epf data appears to be empty. Check that this path is correct and the file isn't empty ", epf_file))
    }
    epf                          <- epf[, -grep("^X...", colnames(epf))]
    epf                          <- t(epf)
    epf                          <- data.frame(epf[, ncol(epf)])
    colnames(epf)                <- "EPF"
    epf$Position_Sample          <- rownames(epf)
    epf[c("Position", "Sample")] <- str_split_fixed(epf$Position_Sample, ": ", 2)
    epf                          <- epf[, c("Position", "Sample", "EPF")]

    # Add experiment name
    exprmnt             <- str_split_fixed(cp_file, assay, 2)[1] # get date and voyage
    exprmnt             <- sub('.*\\/', '', exprmnt)             # get everything after last '/'
    exprmnt             <- gsub('.{1}$', '', exprmnt)            # remove last '_'
    exprmnt_split       <- str_split(exprmnt, "_")
    if (length(exprmnt_split[[1]]) == 3) {
      exprmnt             <- paste0(exprmnt_split[[1]][1], "_", exprmnt_split[[1]][2], exprmnt_split[[1]][3])
    } else if (length(exprmnt_split[[1]]) == 2) {
      exprmnt             <- paste0(exprmnt_split[[1]][1], "_", exprmnt_split[[1]][2])
    } else {
      stop("Your filename isn't valid. The experiment name must contain two or three underscores (e.g., '220912_SWWAV11', or '220912_SWWA_V11').")
    }
    
    epf$Experiment_name <- paste0(exprmnt, "_", assay, "_", curr_plate)
    
    # Import tm and keep only Position, Tm1, and Tm2
    print(paste0("Importing: ", tm_file))
    tm          <- read_delim(tm_file, delim = "\t", skip = 1)
    if (dim(tm)[1] == 0) {
      stop(paste0("tm data appears to be empty. Check that this path is correct and the file isn't empty ", tm_file))
    }
    tm$Position <- tm$Pos
    tm          <- tm[, c("Position", "Tm1", "Tm2")]
    
    # Merge all data, then export results
    output          <- merge(merge(epf, cp, by = "Position"), tm, by = "Position")
    #output_filename <- tail(strsplit(cp_file, "/")[[1]], n = 1)
    #output_filename <- gsub("Cp", "output", output_filename)
    output_filename <- paste0(exprmnt, "_", assay, "_", curr_plate, ".txt")
    print(paste0("Exporting data to: ", output_dir, output_filename))
    
    #Experiment_name	Position	Sample	EPF	Cp	Tm1	Tm2
    output <- output[, c("Experiment_name", "Position", "Sample", "EPF", "Cp", "Tm1", "Tm2")]
    write.table(output, paste0(output_dir, output_filename))
  }
}
