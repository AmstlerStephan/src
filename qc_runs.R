library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(data.table)

runs <- list.dirs(path = "QC/Nanostats", recursive = FALSE)
runs <- runs[ grepl("run", runs)]
print(runs)

str_extract(runs, "run\\d")

run_stats <- list()
stats <- list()
counter <- 1

for (run in runs) {
  for (barcode in list.dirs(paste(run, "stats", sep = "/"), recursive = FALSE)) {
    for (file in list.files(barcode)) {
      nanostat <- read_tsv(paste(barcode, file, sep = "/")) %>%
        transpose(make.names = "Metrics") %>%
        mutate(read = str_split(file, '_', n = 2)[[1]][2],
               barcode = paste('NB', str_sub(basename(barcode), start = -2), sep = ''), 
               run = str_extract(run, "run\\d"))
      
      run_stats[[counter]] <- nanostat
      
      counter <- counter + 1
    }
  }
  
  sample_barcode_overview <-
    fromJSON(paste(str_extract(run, "run\\d"), 'lib/Barcode_Sample_overview.js', sep = "/"))
  
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
