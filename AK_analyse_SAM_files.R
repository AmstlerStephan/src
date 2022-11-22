##########################################################################
### This script is used to read all SAM files created by the           ### 
### make_SAM_files.sh bash script                                      ###
### all unique sequences from the SAM files are calculated and         ###
### saved as a csv_file in the <unique> directory                      ###
### this R script is part of the                                       ###
### collect_unique_sequences.sh script (TO BE DONE)                    ###
##########################################################################

library(plyr)
library(tidyr)
library(readr)
library(tidyverse)


### set result folder (output folder)
args <- commandArgs(TRUE)
if (!is.null(args[1])) {
  input_dir <- args[1]
  print(input_dir)
}

output_dir <- paste(input_dir, "unique", sep = "/")
SAM_dir <- paste(input_dir, "SAM", sep = "/")
SAM_files <- list.files(path = SAM_dir, pattern = "*.sam", full.names = TRUE)

if(!dir.exists(output_dir)) dir.create(output_dir)

trim_seq <- function(CIGAR, sequence) {
  soft_clips <- unlist(str_extract_all(unlist(str_extract_all(CIGAR, "\\d*S")), "\\d+"))
  trimmed_seq <- str_sub(sequence, start = (as.numeric(soft_clips[1]) + 1), end = - (as.numeric(soft_clips[2]) + 1))
}

write_csv_files <- function(i) {
  write.csv(
    unique_sequences[i],
    file = str_replace(str_replace(SAM_files[i], SAM_dir, output_dir), ".sam", ".csv"),
    row.names = FALSE
  )
}

analyze_sam <- function(data) {
  data %>% 
    rowwise() %>%
    mutate(trimmed_seq = trim_seq(X6, X10),
           length_diff = str_length(X10) - str_length(trimmed_seq)) %>% 
    group_by(trimmed_seq) %>% 
    summarize(count = n())
}

unique_sequences <- 
  lapply(lapply(SAM_files,
                read_tsv, col_names = FALSE),
         analyze_sam)

lapply(seq_along(unique_sequences),
       write_csv_files)




