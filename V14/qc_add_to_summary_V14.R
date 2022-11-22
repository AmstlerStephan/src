### check if input values are correct
args <- commandArgs(TRUE)
if(length(args) != 1){
  stop("The script needs exactly 1 argument qc_summary_path)", call.=FALSE)
} 

library(jsonlite)
library(dplyr)
library(tidyverse)
library(janitor)
library(data.table)


### set values for the script
qc_summary_path <- args[1]

runs <- list.files(qc_summary_path, recursive = FALSE, full.names = TRUE, pattern = "run")
summary <- list()
counter <- 1

print(runs)
for (run in runs){
  summary[[counter]] <- read_csv(run, col_names = TRUE)
  
  counter <- counter + 1
  
}

bind_rows(summary) %>% 
  type_convert() %>% 
  write_csv("QC/nanostat_summary_per_run.csv")

            