### check if input values are correct
args <- commandArgs(TRUE)
if(length(args) != 2){
  stop("The script needs exactly 2 arguments (work directory, nanoplot output path, qc output path)", call.=FALSE)
} 

library(jsonlite)
library(dplyr)
library(tidyverse)
library(janitor)
library(data.table)


### set values for the script
run_path <- args[1]
output_path <- args[2]
user_run_id <- basename(run_path)

generell_info <- scan(Sys.glob(paste(run_path, "*.md", sep = "")), what = "character")
json_data <- vector()
key <- generell_info[5]

for(i in 4:length(generell_info)) {
  if(generell_info[i-1] == ","){ key <- generell_info[i] }
  if(generell_info[i-1] == ":"){ value <- generell_info[i] }
  if(i %% 4 == 0 && i != 4){ json_data[[key]] <- value }
  if(generell_info[i] == "}"){ break }
  
}

run_info <- data.frame(as.list(json_data))

### Read Nanoplot run summary

nanoplot <- 
  read_tsv(file = paste( output_path, paste(user_run_id, "txt", sep = "."), sep = "/")) %>% 
  transpose(make.names = "Metrics") %>% 
  mutate( run = user_run_id)

qc <- bind_cols(run_info, nanoplot)

write.csv(qc, paste( paste(output_path, user_run_id, sep = "/"), "csv", sep = "."), row.names = FALSE)
