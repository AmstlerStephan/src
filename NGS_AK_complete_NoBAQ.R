### Script to collect all NGS data, filter it for AK and QiagenLR Polymerase
### Combine it into one dataframe

library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)

NGS_5104 <- read.csv("data_ngs/201803_validation_plasmids_5104_alldata.csv")
NGS_2645 <- read.csv("data_ngs/201804_validation_plasmids_2645_alldata.csv")

### filter for AK samples and select only needed columns
NGS_5104 <- NGS_5104 %>% 
  filter(grepl("AK", as.character(sampleid)) & enzyme == "QiagenLR") %>% 
  select(sample, fragment, gene, pos, ref, variant, major_minor, variant_level)

NGS_2645 <- NGS_2645 %>% 
  filter(grepl("AK", as.character(sampleid)) & enzyme == "QiagenLR") %>% 
  select(sample, fragment, gene, pos, ref, variant, major_minor, variant_level)

NGS_complete <- bind_rows(NGS_2645, NGS_5104)

write.csv(NGS_complete, 'data_ngs/NGS_AK_complete_noBAQ.csv')
write.csv(NGS_5104, 'data_ngs/NGS_AK_5104_noBAQ.csv')
write.csv(NGS_2645, 'data_ngs/NGS_AK_2645_noBAQ.csv')

