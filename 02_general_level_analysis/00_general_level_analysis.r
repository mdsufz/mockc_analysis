## Mock Community analysis

##### Load Packages ####

library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

##### Load table with counts, matches to original lists and False positives ####

all_comm_results_totals <-
  read.csv("01_All_comms_results_R.csv", header = T) #File created on initial_data step

all_comm_results <-
  read.table("02_All_comms_results_R_means.txt",
             header = T,
             sep = "\t") #File created on initial_data step

all_comm_results$Depth <- as.factor(all_comm_results$Depth)

all_comm_results$Abundance_distribution <-
  as.factor(all_comm_results$Abundance_distribution)

##### Evaluate Taxonomic relatedness effect on bin recovery ####

pdf("02_general_level_analysis/bins_per_taxonomy_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3", "4", "5"),
  Taxonomic_distribution %in% c("C", "N", "V"),
  Depth == "60"
) %>% ggplot(.,
             aes(x = Taxonomic_distribution, y = Number_of_bins, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of bins") + xlab("Taxonomic relatedness") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")) + labs(fill = "Evenness distribution")
dev.off()

pdf("02_general_level_analysis/accuracy_per_taxonomy_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3", "4", "5"),
  Taxonomic_distribution %in% c("C", "N", "V"),
  Depth == "60"
) %>% ggplot(.,
             aes(x = Taxonomic_distribution, y = Taxonomies_found, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of original bins found") + xlab("Taxonomic relatedness") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")) + labs(fill = "Evenness distribution")
dev.off()

pdf("02_general_level_analysis/FP_per_taxonomy_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3", "4", "5"),
  Taxonomic_distribution %in% c("C", "N", "V"),
  Depth == "60"
) %>% ggplot(.,
             aes(x = Taxonomic_distribution, y = Number_FP, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of False positive taxonomies") + xlab("Taxonomic relatedness") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")) + labs(fill = "Evenness distribution")
dev.off()

##### Evaluate Sequencing depth on bin recovery ####

pdf("02_general_level_analysis/bins_per_depth_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3"),
  Taxonomic_distribution == "R",
  Depth %in% c("10", "30", "60", "120", "180")
) %>% ggplot(.,
             aes(x = Depth, y = Number_of_bins, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of bins") + xlab("Sequencing depth") +
  scale_fill_manual(values = c("#d7191c", "#fdae61")) + labs(fill = "Evenness distribution")
dev.off()

pdf("02_general_level_analysis/accuracy_per_depth_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3"),
  Taxonomic_distribution == "R",
  Depth %in% c("10", "30", "60", "120", "180")
) %>% ggplot(.,
             aes(x = Depth, y = Taxonomies_found, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of original bins found") + xlab("Sequencing depth") +
  scale_fill_manual(values = c("#d7191c", "#fdae61")) + labs(fill = "Evenness distribution")
dev.off()

pdf("02_general_level_analysis/FP_per_depth_test.pdf")
all_comm_results %>% filter(
  Abundance_distribution %in% c("2", "3"),
  Taxonomic_distribution == "R",
  Depth %in% c("10", "30", "60", "120", "180")
) %>% ggplot(., aes(x = Depth, y = Number_FP, fill = Abundance_distribution)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~ Pipeline, strip.position = "bottom") + theme_classic() + ylab("Average Number of False positive taxonomies") + xlab("Sequencing depth") +
  scale_fill_manual(values = c("#d7191c", "#fdae61")) + labs(fill = "Evenness distribution")
dev.off()

#### Evaluate Evenness profile on bin recovery ####

pdf("02_general_level_analysis/bins_per_evenness_test.pdf")
all_comm_results %>% filter(Abundance_distribution == "1",
                            Taxonomic_distribution == "O",
                            Depth == "60") %>% ggplot(., aes(x = Pipeline, y = Number_of_bins, fill =
                                                               Pipeline)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + ylab("Average Number of bins") + xlab("Pipeline") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9")) + guides(fill =
                                                                            "none")
dev.off()

pdf("02_general_level_analysis/accuracy_per_evenness_test.pdf")
all_comm_results %>% filter(Abundance_distribution == "1",
                            Taxonomic_distribution == "O",
                            Depth == "60") %>% ggplot(., aes(x = Pipeline, y = Taxonomies_found, fill =
                                                               Pipeline)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + ylab("Average Number of original bins found") + xlab("Pipeline") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9")) + guides(fill =
                                                                            "none")
dev.off()

pdf("02_general_level_analysis/FP_per_evenness_test.pdf")
all_comm_results %>% filter(Abundance_distribution == "1",
                            Taxonomic_distribution == "O",
                            Depth == "60") %>% ggplot(., aes(x = Pipeline, y = Number_FP, fill = Pipeline)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() + ylab("Average Number of False positive taxonomies") + xlab("Pipeline") +
  scale_fill_manual(values = c("#d7191c", "#fdae61", "#abd9e9")) + guides(fill =
                                                                            "none")
dev.off()

