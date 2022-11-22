# comparison of NGS variant level vs UMI variant level
# In the NGS dataset AKs are stored as AK_AKXX
# In the mutserve dataset values are stored according to naming of the .bam file
# in this case bcXX_5104_20x.bam

library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)

args <- commandArgs(TRUE)
if (!is.null(args[1])) {
  path <- args[1]
  print(path)
}

NGS <- read.csv('data_ngs/NGS_AK_complete_noBAQ.csv') %>%
  mutate(fragment = as.character(fragment))

barcodes <-
  fromJSON(paste(path, 'lib/Barcode_Sample_overview.js', sep = ''))
UMI <-
  read_tsv(
    paste(
      path,
      '/mutserve/',
      str_sub(path, end = -2),
      '_summary_mutserve.txt',
      sep = ""
    ),
    na = c('', 'NA', '-')
  )

UMI$min_variant <-
  with(UMI, ifelse((`REF` == `MINOR-REV`), `TOP-REV-PERCENT`, `MINOR-REV-PERCENT`))
UMI <- UMI %>% drop_na(`MINOR-REV`) %>%
  separate(SAMPLE,
           c('Barcode', 'Fragment', 'Additional_info', 'Filetype'),
           sep = '[_|.]') %>%
  mutate(
    Barcode = paste('NB', str_sub(Barcode, start = -2), sep = ''),
    Fragment = str_sub(Fragment, start = str_locate(Fragment, '\\d+'))
  ) %>%
  join(barcodes, by = 'Barcode')

UMI_Plasmids <- UMI %>%
  filter(grepl('A_B', Sample)) %>%
  separate(Sample,
           c(NA, NA, 'Percent_A', 'Percent_B', NA),
           sep = '_',
           remove = FALSE) %>%
  mutate(
    Percent_A = as.numeric(Percent_A) / 10 ,
    Percent_B = as.numeric(Percent_B) / 10 ,
    Sample_readable_PL = paste(Percent_A, Percent_B, sep = ':')
  )
UMI_AK <- UMI %>%
  filter(grepl('AK', Sample)) %>%
  mutate(Sample_readable_AK = Sample)

UMI <- bind_rows(UMI_Plasmids, UMI_AK)

UMI$Sample_readable <-
  ifelse(is.na(UMI$Sample_readable_PL),
         UMI$Sample_readable_AK,
         UMI$Sample_readable_PL)


### filter PCR errors

UMI_filtered <- UMI %>% filter(min_variant > 0.005)


### calculate mean and std of deviation to expected value
### especially sd is important to see the consistency of the results
PL_mean_variance <- UMI_filtered %>%
  filter(grepl('A_B', Sample)) %>%
  mutate(type = ifelse(min_variant > 0.5, 'A', 'B')) %>%
  group_by(Sample, type, Percent_A, Percent_B, Sample_readable, Fragment) %>%
  summarise(mean = mean(min_variant),
            SD = sd(min_variant)) %>%
  mutate(divergence = ifelse(type == 'A', Percent_A / 100 - mean, Percent_B / 100 - mean))

NGS_comparison_data <- UMI %>%
  full_join(NGS,
            by = c(
              'Fragment' = 'fragment',
              'Sample' = 'sample',
              'POS' = 'pos'
            )) %>%
  dplyr::rename(
    Position = POS,
    Ref_UMI = REF,
    Minor_UMI = `MINOR-REV`,
    Variant_level_UMI = min_variant,
    Variant_NGS = variant,
    Variant_level_NGS = variant_level,
    Ref_NGS = ref
  ) %>%
  select(
    Sample,
    Fragment,
    Additional_info,
    Position,
    Ref_UMI,
    Minor_UMI,
    Variant_level_UMI,
    Ref_NGS,
    Variant_NGS,
    Variant_level_NGS
  ) %>%
  mutate(
    Variant_level_UMI = coalesce(Variant_level_UMI, 0),
    Variant_level_NGS = coalesce(Variant_level_NGS, 0),
    Variance_level_absolute_difference = Variant_level_NGS - Variant_level_UMI,
    Variance_level_relative_difference = (Variant_level_NGS / Variant_level_UMI - 1)
  ) %>%
  filter(grepl('AK', Sample)) %>%
  filter(Variant_NGS != 'D' &
           Variant_level_NGS < 0.98 & Variant_level_UMI > 0.005) %>%
  filter(is.finite(Variance_level_relative_difference))

### to filter for Mutserve 'errors'
#
# filter(Position < 2471 | Position > 2510)

path <- paste(path, 'results/', sep = '')

write.csv(UMI, paste(path, 'UMI.csv', sep = ''))
write.csv(UMI_filtered, paste(path, 'UMI_filtered.csv', sep = ''))
write.csv(PL_mean_variance, paste(path, 'PL_mean_variance.csv', sep = ''))
write.csv(NGS_comparison_data,
          paste(path, 'NGS_comparison_data.csv', sep = ''))
write.csv(UMI_filtered, paste(path, 'UMI_filtered.csv', sep = ''))

models = NGS_comparison_data %>% group_by(Sample, Fragment) %>% do(model = lm(Variant_level_UMI ~ Variant_level_NGS, data = .))
print(models$model)

##### PLOTS PLOTS PLOTS

if (nrow(NGS_comparison_data) > 2) {
  print(paste("NGS_comparison_data", nrow(NGS_comparison_data)))
  
  test <-
    NGS_comparison_data %>% lm(Variant_level_UMI ~ Variant_level_NGS, data = .)
  print(test)
  print(summary(test)$r.squared)
  # Plotting Variance Levels against each other
  comparison_variant_levels_ngs_umi <- NGS_comparison_data %>%
    #filter(Minor_UMI != Variant_NGS) %>%
    ggplot(aes(x = Variant_level_UMI, y = Variant_level_NGS)) +
    facet_wrap(vars(Sample, Fragment)) +
    geom_smooth(method = 'lm', show.legend = TRUE) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    geom_point(position = 'jitter') +
    labs(x = 'relative variance level UMI',
         y = 'relative variance level NGS',
         title = 'Variant levels of both Sequencing technology') +
    annotate(
      'text',
      x = 0.8,
      y = 0.05,
      label = paste0('R Squared = ', round(summary(test)$r.squared, digits = 4))
    ) +
	xlim(0,1)
  
  
  
  # Compare variance levels of different Sequencing methods
  comparison_variant_levels_ngs_umi_per_position <-
    NGS_comparison_data %>%
    #filter(Variant_level_UMI != 0 && Variant_level_NGS != 0) %>%
    ggplot(aes(x = Position)) +
    facet_wrap(vars(Sample, Fragment)) +
    geom_point(aes(y = Variant_level_NGS, color = 'NGS'),
               show.legend = TRUE) +
    geom_point(aes(y = Variant_level_UMI, color = 'UMI'),
               show.legend = TRUE) +
    labs(x = 'position relative to reference sequence',
         y = 'relative variant Level',
         title = 'Variant levels per position and Sequencing technology')
  
  # Plot showing mean difference of variance level of both Sequencing Methods
  absolute_difference_variance_level <- NGS_comparison_data %>%
    ggplot(aes(x = Position, y = Variance_level_absolute_difference)) +
    facet_wrap(vars(Sample, Fragment)) +
    geom_point() +
    geom_label(
      data = NGS_comparison_data %>%
        filter(
          Variance_level_absolute_difference > 0.055 |
            Variance_level_absolute_difference < -0.055
        ),
      # Filter data first
      aes(label = Position)
    ) +
    labs(x = 'Position',
         y = 'Absolute Differences in Variant level [ NGS - UMI ]',
         title = 'Absolute Difference in Variant Levels per Position')
  
  # Plot showing relative difference of variance level of both Sequencing Methods
  relative_difference_variance_level <-  NGS_comparison_data %>%
    ggplot(aes(x = Position, y = Variance_level_relative_difference)) +
    facet_wrap(vars(Sample, Fragment)) +
    geom_point() +
    geom_label(
      data = NGS_comparison_data %>%
        filter(
          Variance_level_relative_difference > 1 |
            Variance_level_relative_difference < -1
        ),
      # Filter data first
      aes(label = Position)
    ) +
    labs(x = 'Position',
         y = 'Relative Differences in Variant level [ NGS / UMI - 1]',
         title = 'Relative Difference in Variant Levels per Position')
  
  # Comparing relative and absolute difference in Variant level
  absolute_vs_relative_difference_variance_level <-
    NGS_comparison_data %>%
    ggplot(
      aes(x = Variance_level_absolute_difference, y = Variance_level_relative_difference)
    ) +
    facet_wrap(vars(Sample, Fragment)) +
    geom_point() +
    labs(x = 'Absolute Differences in Variant level [ NGS - UMI ]',
         y = 'Relative Differences in Variant level [ NGS / UMI - 1]',
         title = 'Comparison of absolute and relative Differences in Variant levels')
  
  
  ggsave(
    'NGS_comparison_variant_levels_ngs_umi.jpg',
    path = path,
    comparison_variant_levels_ngs_umi,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  ggsave(
    'NGS_comparison_variant_levels_ngs_umi_per_position.jpg',
    path = path,
    comparison_variant_levels_ngs_umi_per_position,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  ggsave(
    'NGS_absolute_difference_variance_level.jpg',
    path = path,
    absolute_difference_variance_level,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  ggsave(
    'NGS_relative_difference_variance_level.jpg',
    path = path,
    relative_difference_variance_level,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  ggsave(
    'NGS_absolute_vs_relative_difference_variance_level.jpg',
    path = path,
    absolute_vs_relative_difference_variance_level,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  
}

##### PLOTS PLASMIDS

PL_mean_variance <- PL_mean_variance %>% filter(!grepl('AK', Sample) & Sample != 'A_B_950_50_2645' & Sample != 'A_B_990_10_2645')

if (nrow(PL_mean_variance) > 2) {
  print(paste("PL_mean_variance:", nrow(PL_mean_variance)))
  
  plot_PL_mean_variance <- ggplot(PL_mean_variance,
                                  aes(x = Sample_readable, y = divergence)) +
    facet_wrap(vars(type, Fragment)) +
    geom_bar(stat = 'identity') +
    labs(x = 'Sample acronym',
         y = 'Mean divergence compared to expected value of mixture',
         title = 'Mean divergence of variant levels') +
    theme(axis.text.x = element_text(angle = 8, vjust = 0.5))
  
  ggsave(
    'PL_mean_variance.jpg',
    plot_PL_mean_variance,
    path = path,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  
}

if (nrow(UMI_filtered %>% filter(grepl('A_B', Sample))) > 2) {
  print(paste("UMI filtered_A_B:", nrow(UMI_filtered %>% filter(
    grepl('A_B', Sample)
  ))))
  
  plot_levels <- ggplot(
    UMI_filtered %>%
      filter(grepl('A_B', Sample)) %>%
      mutate(type = ifelse(min_variant > 0.5, 'A', 'B')),
    aes(y = min_variant, fill = type)
  ) +
    facet_wrap(vars(Sample_readable, Fragment)) +
    geom_boxplot()
  
  ggsave(
    'PL_variance_levels_boxplot.jpg',
    plot_levels,
    path = path,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  
}

UMI_filtered <- UMI_filtered %>% filter(!grepl('AK', Sample) & Sample != 'A_B_950_50_2645' & Sample != 'A_B_990_10_2645')

if (nrow(UMI_filtered) > 2) {
  print(paste("UMI filtered:", nrow(UMI_filtered)))
  
  plot_variance_level_per_sample <- ggplot(UMI_filtered,
                                           aes(x = POS, y = min_variant)) +
    facet_wrap(vars(Sample_readable, Fragment)) +
    geom_point() +
    labs(x = 'Position compared to reference genome',
         y = 'Variant level',
         title = 'Variant level across the whole Amplicon')
  
  plot_variance_level_per_sample_zoomed <- ggplot(UMI_filtered,
                                                  aes(x = POS, y = min_variant)) +
    facet_wrap(vars(Sample_readable, Fragment)) +
    geom_point() +
    coord_cartesian(ylim = c(0.00, 0.1)) +
    labs(x = 'Position compared to reference genome',
         y = 'Variant level',
         title = 'Variant level across the whole amplicon')
  
  ggsave(
    'ALL_variance_level_per_sample_zoomed.jpg',
    plot_variance_level_per_sample_zoomed,
    path = path,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  ggsave(
    'ALL_variance_level_per_sample.jpg',
    plot_variance_level_per_sample,
    path = path,
    width = 11.35,
    height = 8.42,
    dpi = 300
  )
  
}
