library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
library(data.table)

barcodes_path <-
  list.dirs("run7/Nanostats/stats", recursive = FALSE)
stats <- list()
counter <- 1

for (i in barcodes_path) {
  for (j in list.files(i)) {
    
    nanostat <- read_tsv(paste(i, j, sep = "/")) %>% 
      transpose(make.names = "Metrics") %>%
      mutate(read = str_split(j, '_', n = 2)[[1]][2],
             barcode = basename(i))

    stats[[counter]] <- nanostat
    
    counter <- counter + 1
  }
}

stats <- bind_rows(stats) %>% 
  type_convert() %>% 
  group_by(barcode) %>%
  filter(number_of_reads > 10000) %>% 
  mutate(max = max(mean_qual),
         min = min(mean_qual),
         diff = max(mean_qual) - min(mean_qual))

stats %>% 
  ggplot(aes(barcode, mean_qual, color = read)) +
  geom_point() +
  geom_segment(aes(x=barcode, xend=barcode, y=min, yend=max)) +
  coord_cartesian(ylim= c(10,13)) +
  scale_y_continuous(breaks = seq(10,13, 1))
  
