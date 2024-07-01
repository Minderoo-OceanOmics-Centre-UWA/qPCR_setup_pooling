##########################################################
# Libraries
##########################################################

library(readr)
library(stringr)


##########################################################
# Variables
##########################################################

input_dir  <- "test_data/input/cp_epf_tm_files/"
output_dir <- "test_data/output/"
assays     <- c("16S", "MiFish")


##########################################################
# Main
##########################################################

for (assay in assays) {
  cp_files    <- Sys.glob(paste0(input_dir, "*", assay, "*Cp.txt"))
  plate_count <- length(cp_files)
  
  for (plate_num in 1:plate_count) {
    curr_plate <- paste0("Plate", plate_num)
    print(paste0(assay, ": ", curr_plate))
    
    cp_file  <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*Cp.txt"))
    epf_file <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*EPF.txt"))
    tm_file  <- Sys.glob(paste0(input_dir, "*", assay, "_", curr_plate, "*Tm.txt"))
    
    # Import cp and keep only Position, and Cp
    print(paste0("Importing: ", cp_file))
    cp          <- read_delim(cp_file, delim = "\t", skip = 1)
    cp$Position <- cp$Pos
    cp          <- cp[, c("Position", "Cp")]
    
    # Import epf, remove weird columns, transpose, keep final EPF plus Position and Sample
    print(paste0("Importing: ", epf_file))
    epf                          <- read_delim(epf_file, delim = "\t")
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
    epf$Experiment_name <- paste0(exprmnt, "_", assay, "_", curr_plate)
    
    # Import tm and keep only Position, Tm1, and Tm2
    print(paste0("Importing: ", tm_file))
    tm          <- read_delim(tm_file, delim = "\t", skip = 1)
    tm$Position <- tm$Pos
    tm          <- tm[, c("Position", "Tm1", "Tm2")]
    
    # Merge all data, then export results
    output          <- merge(merge(epf, cp, by = "Position"), tm, by = "Position")
    output_filename <- tail(strsplit(cp_file, "/")[[1]], n = 1)
    output_filename <- gsub("Cp", "output", output_filename)
    print(paste0("Exporting data to: ", output_dir, output_filename))
    write.table(output, paste0(output_dir, output_filename))
  }
}
