library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)


### define parameters 
STR_start <- 2472
STR_end <- 2505
UMI_cutoff <- 0.0085

args <- commandArgs(TRUE)
if (!is.null(args[1])) {
  path <- args[1]
  print(path)
}

run <- str_extract(path, "run\\d")

mutserve_summary <-
  read_tsv(
    paste(
      path,
      '/mutserve/',
      run,
      '_summary_mutserve.txt',
      sep = ""
    ),
    na = c('', 'NA', '-')
  )

barcodes <-
  read.csv("QC/nanostats_barcode_all_runs.csv") %>% 
  filter(read == "raw") %>% 
  select(barcode:max, number_of_reads) %>% 
  mutate(Q_score = max)

run_info <- 
  read_csv("info/Status_overview.csv") %>% 
  select(Acronym, Run_directory, fastq_pass_reads, Fragment)

NGS <- read_csv("data_ngs/NGS_AK_complete_noBAQ.csv") %>% 
  filter(pos < STR_start | pos > STR_end)

### filter mutserve data 
### filter for full conversions (called variant is not the reference AND has no minor variant level OR minor variant level is below a certain threshold)
mutserve_raw_full_conversions <- mutserve_summary %>% 
  filter(REF != `TOP-REV` & (is.na(`MINOR-REV`) | `MINOR-REV-PERCENT` < UMI_cutoff)) %>% 
  mutate(Variant_level_UMI = 1, 
         Variant_UMI = `TOP-REV`)

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
mutserve_raw_variants <- mutserve_summary %>% 
  drop_na(`MINOR-REV`) %>% 
  mutate(Variant_level_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV-PERCENT`, `MINOR-REV-PERCENT`),
         Variant_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV`, `MINOR-REV`))

mutserve_combined <- bind_rows(mutserve_raw_full_conversions, mutserve_raw_variants)

### Join Barcodes and mutserve data

UMI <- mutserve_combined %>%
  separate(SAMPLE,
           c('Barcode', 'Fragment', 'Additional_info', 'Filetype'),
           sep = '[_|.]') %>%
  mutate(
    barcode = paste('NB', str_sub(Barcode, start = -2), sep = ''),
    Fragment = as.numeric(str_sub(Fragment, start = str_locate(Fragment, '\\d+'))),
    run = run
  ) %>%
  inner_join(barcodes, by = c('barcode', 'run')) %>% 
  dplyr::rename(pos = POS, 
                sample = Sample,
                fragment = Fragment)

UMI_Plasmids <- UMI %>%
  filter(grepl('A_B', sample)) %>%
  separate(sample,
           c(NA, NA, 'Percent_A', 'Percent_B', NA),
           sep = '_',
           remove = FALSE) %>%
  mutate(
    Percent_A = as.numeric(Percent_A) / 10 ,
    Percent_B = as.numeric(Percent_B) / 10 ,
    Sample_readable_PL = paste(Percent_A, Percent_B, sep = ':')
  )
UMI_AK <- UMI %>%
  filter(grepl('AK', sample) | grepl('5400', sample)) %>%
  mutate(Sample_readable_AK = sample)

UMI <- bind_rows(UMI_Plasmids, UMI_AK)

UMI$Sample_readable <-
  ifelse(is.na(UMI$Sample_readable_PL),
         UMI$Sample_readable_AK,
         UMI$Sample_readable_PL)


UMI_filtered <- UMI %>% 
  filter(Variant_level_UMI > UMI_cutoff) %>%
  filter(pos < STR_start | pos > STR_end)

### getting a dataframge with AKs according fragments included in the sequencing run 
UMI_groups_AK <- UMI_AK %>% group_by(sample, fragment) %>% summarize()

### the groups are used to filter the NGS data before joining both dataframes
### exclude all samples that are not covered with the UMI run

NGS <- NGS %>%
  filter(sample %in% UMI_groups_AK$sample & fragment %in% UMI_groups_AK$fragment)

NGS <- NGS %>% 
  full_join(UMI_AK,
            by = c('fragment', 'sample', 'pos')) %>%
  dplyr::rename(
    Position = pos,
    Ref_UMI = REF,
    Variant_NGS = variant,
    Variant_level_NGS = variant_level,
    Ref_NGS = ref, 
    num_of_consensus_sequences = `COV-TOTAL`
  ) %>%
  select(
    sample,
    fragment,
    Additional_info,
    Position,
    Ref_UMI,
    Variant_UMI,
    Variant_level_UMI,
    Ref_NGS,
    Variant_NGS,
    Variant_level_NGS, 
    run, 
    number_of_reads,
    num_of_consensus_sequences, 
    Q_score
  ) %>%
  mutate(
    Variant_level_UMI = coalesce(Variant_level_UMI, 0),
    Variant_level_NGS = coalesce(Variant_level_NGS, 0),
    Variance_level_absolute_difference = Variant_level_NGS - Variant_level_UMI,
    # Variance_level_relative_difference = (Variant_level_NGS / Variant_level_UMI - 1)
  )

NGS$run <- run

NGS <- NGS %>% 
  inner_join(run_info, by = c('sample' = 'Acronym', 'run' = 'Run_directory', 'fragment' = 'Fragment')) 

NGS_filtered <- NGS %>% 
  filter(!is.na(Variant_NGS) | (Variant_level_UMI > UMI_cutoff | Variant_level_UMI == 0)) %>%
  filter(Position < STR_start | Position > STR_end) %>%
  filter(Variant_NGS != 'D')

path <- paste(path, 'results/', sep = '')

write.csv(UMI, paste(path, 'UMI.csv', sep = ''))
write.csv(UMI_filtered, paste(path, 'UMI_filtered.csv', sep = ''))
write.csv(mutserve_combined, paste(path, 'mutserve_combined.csv', sep = ''))
write.csv(NGS, paste(path, 'NGS.csv', sep = ''))
write.csv(NGS_filtered, paste(path, 'NGS_filtered.csv', sep = ''))




