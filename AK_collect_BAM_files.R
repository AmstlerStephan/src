##########################################################################
### This script collects all AK_samples of one run and combines them   ###
### into a standardized folder structure (out_dir/BAM/<SAMPLE>/)       ###
### files are collected from the run*/ont_pl/barcode*/align/ folder    ###
### AK_samples plus according barcodes are found in the                ###
### info/Barcode_Sample_overview_summary.js file , created with        ###
### INSERT.sh script                                                   ###
### only the newest run of an AK is used for analysis                  ###
##########################################################################

library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)

### set result folder (output folder)
args <- commandArgs(TRUE)
if (!is.null(args[1])) {
  out_dir <- args[1]
  print(out_dir)
}

out_dir <- paste(out_dir, "BAM", sep = "/")
bam_type <- "final"

json <- read_json("info/Barcode_Sample_overview_summary.js", simplifyVector = TRUE)

AK_extract_barcodes <- function() {
  list <- list()
  counter <- 1
  for (run in json$Runs) {
    list[[counter]] <- run %>% drop_na() %>%
      # filter(grepl("AK", Sample)) %>%
      mutate(run = names(json$Runs)[counter])
    counter = counter + 1
  }
  return(list)
}

collect_AK_samples <- function(data) {
  for (row in 1:nrow(data)) {
    barcode <- data[row, 1]
    sample <- data[row, 2]
    run <- data[row, 3]
    
    if (!dir.exists(run)) {
      next
    }
    
    target_dir <- paste(out_dir, run, sample, sep = "/")
    src_dir <- paste(run, "ont_pl", barcode, "align", sep = "/")
    files_to_copy <- list.files(path = src_dir, pattern = bam_type, full.names = TRUE)
    
    print(src_dir)
    print(files_to_copy)

    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    
    print(target_dir)
    
    file.copy(files_to_copy, target_dir, overwrite = TRUE)
  }
}


AK_barcodes <- bind_rows(AK_extract_barcodes()) 
collect_AK_samples(AK_barcodes)

