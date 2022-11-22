library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)

sample_sets <- c("SAPHIR")
SAPHIR_samples <- c("4612|4901|4451|4624|4864|5248|5538|5400")
AK_samples <- c("AK03|AK07|AK14|AK17|AK33")
enzymes <- c("QiagenLR")

NGS_SAPHIR <- read_csv("data_ngs/20221122_201803_pop_variants_noBAQ.csv")
NGS_5104 <- read.csv("data_ngs/20221122_201803_validation_plasmids_5104_alldata.csv")
NGS_2645 <- read.csv("data_ngs/20221122_201804_validation_plasmids_2645_alldata.csv")
mutation_classification <- read_csv("data_ngs/20221122_plasmid_expected_muts.csv") # TOD update mutations with subtype A and Type C 

### filter for AK samples and select only needed columns
NGS_5104_filtered <- NGS_5104 %>% 
  filter(grepl(AK_samples , sample) & grepl(enzymes, enzyme)) %>% 
  select(sample, fragment, gene, pos, ref, variant, major_minor, variant_level)

NGS_2645_filtered <- NGS_2645 %>% 
  filter(grepl(AK_samples, sample) & grepl(enzymes, enzyme)) %>% 
  select(sample, fragment, gene, pos, ref, variant, major_minor, variant_level)

NGS_SAPHIR_filtered <- NGS_SAPHIR %>% 
  filter( grepl(sample_sets, sample_sets)) %>% 
  filter( grepl(SAPHIR_samples, id)) %>%
  mutate( 
    sample = id,
    fragment = as.integer(str_sub(runid, start = -6, end = -3)),
    gene = "LPA",
    variant_level = as.numeric(variant_level)
          ) %>% 
  select(sample, fragment, gene, pos, ref, variant, major_minor, variant_level)

NGS_complete <- bind_rows(NGS_2645_filtered, NGS_5104_filtered, NGS_SAPHIR_filtered)

NGS_complete_mutation_classification <- NGS_complete %>% 
  left_join(mutation_classification, by = c("pos" = "Position" , "fragment"))

write_csv(NGS_complete, file = "data_ngs/20221122_NGS_reference_data_SAPHIR.csv")
write_csv(NGS_complete_mutation_classification, file = "data_ngs/20221122_NGS_reference_data_SAPHIR_mutation_classification.csv")
