library(plyr)
library(tidyr)
library(readr)
library(tidyverse)
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
library(Rsamtools)
library(GenomicAlignments)

input_dir <- "subsampling/"
mutserve_dir <- paste(input_dir, "mutserve", sep = "")
ont_pl_dir <- paste(input_dir, "ont_pl", sep = "")
figure_dir <- "figures"

get_sample_name <- function(x, pos){
    x[pos]
}
get_named_sample_list <- function(dir) {
  samples <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
  print(samples)
  samples_names <- lapply(str_split(samples, "/"), FUN = get_sample_name, 3)
  samples_list <- as.list(samples)
  names(samples_list) <- samples_names
  return(samples_list)
}
read_tsv_mutserve <- function(name, sample) {
  read_tsv(sample[name]) %>% 
    mutate(subsample = name)
}
merge_tsvs_mutserve <- function(sample) {
  bind_rows(lapply(names(sample), read_tsv_mutserve, sample)) %>% na_if("-")
}
get_nested_folder_structure <- function(x, pattern, pos) {
  subsample <- list.dirs(x, full.names = TRUE, recursive = TRUE)
  files <- list.files(subsample, pattern = pattern, full.names = TRUE)
  subsample_names <- sapply(str_split(files, "/"), FUN = get_sample_name, pos, simplify = TRUE)
  files <- as.list(files)
  names(files) <- subsample_names
  return(files)
}
### filter mutserve data 
### filter for full conversions (called variant is not the reference AND has no minor variant level)
get_mutserve_raw_full_conversions <- function(data){ 
  data %>% 
    filter(REF != `TOP-REV` & is.na(`MINOR-REV`))  %>% 
    mutate(Variant_level_UMI = 1, 
           Variant_UMI = `TOP-REV`)
}

### drop all NA values, where no minor variant was found (Either full conversion or no variant at that position)
### Variant level can be over 50% -> take Percentage and variant accordingly
get_mutserve_raw_variants <- function(data){ 
  data %>% 
    drop_na(`MINOR-REV`) %>% 
    mutate(Variant_level_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV-PERCENT`, `MINOR-REV-PERCENT`),
           Variant_UMI = ifelse((`REF` == `MINOR-REV`), `TOP-REV`, `MINOR-REV`))
}

get_bam_reads_summary_statistics <- function(ont_pl_directory_structure_list) {
  mapply(ont_pl_directory_structure_list, names(ont_pl_directory_structure_list), FUN = function(sample, sample_name) {
    bind_rows(mapply(sample, names(sample), FUN = function(subsample, subsample_name) {
      summary <- tibble(sample = sample_name,
                        subsample = subsample_name)
      add_column(summary, countBam(subsample))
    }, SIMPLIFY = FALSE
    ))
  }, SIMPLIFY = FALSE
  )
} 


mutserve_directory_structure_list <- lapply(
  get_named_sample_list(mutserve_dir), get_nested_folder_structure, "raw.txt", 4
)

mutserve_raw_summary <- lapply(mutserve_directory_structure_list, merge_tsvs_mutserve)

summary_statistics_mutserve <- tibble(sample = character(), 
                             subsample = character(), 
                             num_of_mutations = numeric(),
                             mean_mut_level = numeric(), 
                             median_mut_level = numeric())

create_plots_mutserve <- lapply(seq_along(mutserve_raw_summary), function(i, samples) {
  name <- names(samples)[i]
  sample <- samples[[i]]
  
  mutserve_combined <- bind_rows(get_mutserve_raw_full_conversions(sample), get_mutserve_raw_variants(sample))
  
  variant_level_per_pos_per_subsample <- mutserve_combined %>% filter(Variant_level_UMI > 0.0085) %>% ggplot(aes(POS, Variant_level_UMI)) +
    facet_wrap(vars(subsample)) +
    geom_point() +
    labs(
      x = 'position relative to reference sequence',
      y = 'relative variant Level',
      title = 
        'Variant levels per position and Sequencing technology'
    ) +
    coord_cartesian(ylim = c(0, 0.1)) +
    scale_y_continuous(breaks = seq(0, 0.1, by = 0.02))
  
  
  variant_level_per_position <- mutserve_combined %>% filter(Variant_level_UMI > 0.0085) %>%  ggplot(aes(POS, Variant_level_UMI, color = subsample)) +
    geom_point() +
    labs(
      x = 'position relative to reference sequence',
      y = 'relative variant Level',
      title = 'Variant levels per position and Sequencing technology'
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1))
  
  
  out_dir <- paste(mutserve_dir, figure_dir, name, sep = "/")
  dir.create(out_dir, recursive = TRUE)
  
  ggsave(filename = "variant_level_per_pos_per_subsample", plot = variant_level_per_pos_per_subsample, path = out_dir, device = "jpeg")
  ggsave(filename = "variant_level_per_position", plot = variant_level_per_position, path = out_dir, device = "jpeg")
}, mutserve_raw_summary)

test_sample <- mutserve_raw_summary$A_B_950_50_2645
test_mutserve_combined <- bind_rows(get_mutserve_raw_full_conversions(test_sample), get_mutserve_raw_variants(test_sample))


ont_pl_directory_structure_list <- lapply(
  get_named_sample_list(ont_pl_dir), get_nested_folder_structure, "final.bam$", 4
)


out_dir <- paste(ont_pl_dir, figure_dir, sep = "/")
dir.create(out_dir, recursive = TRUE)

write_tsv(bind_rows(get_bam_reads_summary_statistics_mapply(ont_pl_directory_structure_list)
), file = paste(out_dir, "final_bam_file_summary_statistics.tsv", sep = "/")) 
