### check if input values are correct
args <- commandArgs(TRUE)
if(length(args) != 2){
  stop("The script needs exactly 2 arguments (work directory, nanoplot output path, qc output path)", call.=FALSE)
} 

#install.packages("jsonlite")
#install.packages("dplyr")
#install.packages("tidyverse")
#install.packages("janitor")


library(jsonlite)
library(dplyr)
library(tidyverse)
library(janitor)


### set values for the script
run_path <- args[1]
output_path <- args[2]
nanoplot_path <- paste(output_path, "/Nanoplot", sep = "" )
user_run_id <- basename(run_path)

#print(run_path)
#print(nanoplot_path)
#print(output_path)
#print(user_run_id)

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

nanoplot <- read_tsv(paste(nanoplot_path, "NanoStats.txt", sep = "/"))
nanoplot <- t(nanoplot) 
nanoplot <- nanoplot %>% row_to_names(row_number = 1)


qc <- tibble(user_run_id)
qc <- bind_cols(qc, run_info, nanoplot)

write.csv(qc, paste(output_path, paste(user_run_id, "QC", sep = "_"), sep = "/"), row.names = 1)
