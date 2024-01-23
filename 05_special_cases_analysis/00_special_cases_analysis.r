#0.2) Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)

#0.3) Load data

full_results_v2 <-
  read.csv("joao_scripts/missingfiles/full_results_v2.csv", header = T)

total_coverages <- read.csv("joao_scripts/missingfiles/coverage_table_rstudio.csv",header=T,sep=",")

####
#### Analysis of special cases ####
####

special_cases_list_A <-
  c("Brucella melitensis",
    "Rhodobacter sphaeroides",
    "Vibrio cholerae")

special_cases_list_B <-
  c("Brucella ovis", "Paracoccus denitrificans", "Vibrio fluvialis")

special_cases_list_C <-
  c("Bacillus subtilis",
    "Borrelia burgdorferi",
    "Streptomyces griseus")

special_cases_all_comms <-
  c(
    "Brucella melitensis",
    "Rhodobacter sphaeroides",
    "Vibrio cholerae",
    "Brucella ovis",
    "Paracoccus denitrificans",
    "Vibrio fluvialis",
    "Bacillus subtilis",
    "Borrelia burgdorferi",
    "Streptomyces griseus"
  )


special_cases <- full_results_v2

special_cases %>%
  arrange(Pipeline) %>%
  select(
    Community,
    Evenness.distribution,
    Taxonomic.distribution,
    Sequencing.Depth..million.reads.,
    Species
  )

complete_cases <- special_cases %>%
  group_by(Pipeline,
           Taxonomic.distribution,
           Evenness.distribution,
           Sequencing.Depth..million.reads.) %>%
  distinct(
    Evenness.distribution,
    Taxonomic.distribution,
    Sequencing.Depth..million.reads.,
    Pipeline,
    Species,
    .keep_all = F
  ) %>%
  ungroup()


#### Analysis of TP low coverage and FP high coverage ####

complete_cases <- as.data.frame(complete_cases)
evenness_list <- unique(complete_cases$Evenness.distribution)
tax_list <- unique(complete_cases$Taxonomic.distribution)
sequencing_list <- unique(complete_cases$Sequencing.Depth..million.reads.)
pipeline_list <- unique(complete_cases$Pipeline)

list_results <- list()

for (i in 1:length(evenness_list)) {
  for (j in 1:length(tax_list)) {
    for (k in 1:length(sequencing_list)) {
      for (l in 1:length(pipeline_list)) {
        items_names <-
          paste(
            toString(evenness_list[i]),
            toString(tax_list[j]),
            toString(sequencing_list[k]),
            toString(pipeline_list[l]),
            sep = "_"
          )
        temp_list <-
          complete_cases[complete_cases$Evenness.distribution == evenness_list[i] &
                           complete_cases$Taxonomic.distribution == tax_list[j] &
                           complete_cases$Pipeline == pipeline_list[l] &
                           complete_cases$Sequencing.Depth..million.reads. == sequencing_list[k], ]
        temp_species <- temp_list$Species
        special_hits_list <-
          list(temp_species[temp_species %in% special_cases_all_comms])
        names(special_hits_list) <- items_names
        list_results <- append(list_results, special_hits_list)
      }
    }
  }
}

list_results_special_cases <-
  list_results[lapply(list_results, length) > 0] # Remove empty lists


final_matrix <- matrix(ncol = 2)

for (l in 1:length(list_results_special_cases)) {
  row_numbers <- list_results_special_cases[[l]]
  temp_matrix <- matrix(nrow = length(row_numbers), ncol = 2)
  temp_matrix[, 2] <- matrix(unlist(list_results_special_cases[l]))
  temp_matrix[, 1] <- names(list_results_special_cases[l])
  final_matrix <- rbind(final_matrix, temp_matrix)
}
final_matrix <- final_matrix[-1, ]

write.table(final_matrix, file = "special_cases_recovered.tsv", sep = "\t")

special_cases_counts <- as.matrix(table(final_matrix[, 1]))

write.table(special_cases_counts, file = "total_number_of_special_cases_per_comm.tsv", sep =
              "\t")

total_coverages %>% group_by(
  evenness,
  taxonomy,
  depth,
  pipeline
) %>%
  summarise(Ave_cov = mean(coverage),
            ave_rel = mean(rel_abd))


cov_analysis <- total_coverages %>%
  group_by(eve_class, taxonomy, depth, pipeline) %>%
  summarise(
    n = n(),
    mean_cov = mean(coverage),
    sd_cov = sd(coverage),
    mean_rel = mean(rel_abd),
    sd_rel = sd(rel_abd),
  ) %>%
  ungroup() %>%
  print(n = Inf)

write.table(cov_analysis, "all_mean_coverages.csv")