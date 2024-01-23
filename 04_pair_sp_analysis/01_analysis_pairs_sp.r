#) Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)

#### Analysis of pairs of species  ####

species_pairs <-
  read.csv(
    "04_pair_sp_analysis/05_pair_recov_new_triplets_summary.csv",
    header = T,
    sep = ","
  ) # Generated on 00_preprocess_pairs_sp.r

species_pairs %>%
  arrange(Pipeline) %>%
  select(Taxonomic_distribution,
         Evenness_distribution,
         Found.percentage)

mean_pairs <- species_pairs %>%
  group_by(Pipeline, Taxonomic_distribution, Evenness_distribution) %>%
  summarise(
    n = n(),
    mean = mean(Found.percentage),
    sd = sd(Found.percentage)
  ) %>%
  ungroup() %>%
  print(n = Inf)

#write.table(mean_pairs, file = "04_pair_sp_analysis/mean_percentages_species_pairs.tsv", sep =
#              "\t")

# do stat test

stat.test <- species_pairs %>%
  group_by(Evenness_distribution, Taxonomic_distribution) %>%
  rstatix::t_test(Found.percentage ~ Pipeline, p.adjust.method = "bonferroni")

# Remove unnecessary columns and display the outputs
taxonomy_test <- stat.test %>% select(-.y.,-statistic,-df) %>%
  print(n = Inf)

write.table(taxonomy_test, file = "04_pair_sp_analysis/Stat_test_pairs_taxonomy_statistics_results_pipelines.csv")

pdf("04_pair_sp_analysis/species_pairs_pipelines.pdf") 

myplot <- ggpubr::ggboxplot(
  species_pairs,
  x = "Pipeline",
  y = "Found.percentage",
  fill = "Pipeline",
  palette = "npg",
  legend = "none"
) +
  facet_wrap( ~ Evenness_distribution * Taxonomic_distribution)
# Add statistical test p-values
stat.test <- stat.test %>% rstatix::add_xy_position(x = "Pipeline")
myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p")

dev.off()
