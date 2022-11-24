library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(data.table)

runs <- list.dirs(path = "QC/Nanostats/", recursive = FALSE)
runs <- runs[ grepl("run*", runs)]
print(runs)

print(str_extract(runs, "run\\d*_*[[:alnum:]]*"))

run_stats <- list()
stats <- list()
counter <- 1

for (run in runs) {
  stats_dir <- paste(run, "stats", sep = "/")
  for (barcode_dir in list.dirs(stats_dir, recursive = FALSE)) {
    for (file in list.files(barcode_dir, full.names = TRUE)) {
      file_name <- basename(file)
      read_type <- str_split(file_name, "_", simplify = TRUE)[2]
      barcode <- basename(barcode_dir)
      nanostat <- read_tsv(file) %>%
        transpose(make.names = "Metrics") %>%
        mutate(read = read_type,
               barcode = paste('NB', str_sub(barcode, start = -2), sep = ''), 
               run = str_extract(run, "run\\d*_*[[:alnum:]]*"))
  
      run_stats[[counter]] <- nanostat
      
      counter <- counter + 1
    }
  }
  
  sample_barcode_overview <-
    fromJSON(paste(str_extract(run, "run\\d*_*[[:alnum:]]*"), 'lib/Barcode_Sample_overview.js', sep = "/"))
  
  run_stats <- bind_rows(run_stats) %>% 
    inner_join(sample_barcode_overview, by = c("barcode" = "Barcode"))
  
  if( is_empty(stats)){
    stats <- bind_rows(run_stats)
  } else {
    stats <- bind_rows(stats, run_stats)
  }
  
  counter <- 1
  run_stats <- list()
  
}

stats <- bind_rows(stats) %>% 
  type_convert() %>% 
  group_by(barcode) %>%
  filter(number_of_reads > 10000) %>% 
  mutate(max = max(mean_qual),
         min = min(mean_qual),
         diff = max(mean_qual) - min(mean_qual))

write.csv(stats, "QC/nanostats_barcode_all_runs.csv", row.names = FALSE)


compare_qual_after_filtering_lollipop <- stats %>% 
  ggplot(aes(Sample, mean_qual, color = read)) +
  facet_wrap(vars(run)) +
  geom_point() +
  geom_segment(aes(x=Sample, xend=Sample, y=min, yend=max)) +
  coord_cartesian(ylim= c(10,20)) +
  scale_y_continuous(breaks = seq(10,20, 1))

ggsave(
  "compare_qual_after_filtering_lollipop.jpg",
  compare_qual_after_filtering_lollipop,
  path = "QC/"
)
