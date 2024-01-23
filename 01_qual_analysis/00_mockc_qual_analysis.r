## Mock Community analysis

##### Load Packages ####

library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

#### Load Data ####
#File created after the run from 00_initial_data
results_all_bin_counts <-
  read.csv("01_All_comms_results_R.csv",
           header = T,
           sep = ",")

##### Evaluate Abundance distribution ####

AD_analysis <-
  results_all_bin_counts[results_all_bin_counts$Depth == 60 &
                           results_all_bin_counts$Taxonomic_distribution == "O", ]
AD_results <-
  AD_analysis %>% group_by(Pipeline) %>% summarise(
    Average_bins = mean(Number_of_bins),
    Percentage_taxa_found = (mean(Taxonomies_found) / 42 * 100),
    Average_FP = mean(Number_FP)
  )

#write.table(AD_results,file="AD_results.tsv",sep = "\t")

##### Evaluate Depth influence ####

Depth_analysis <-
  results_all_bin_counts[results_all_bin_counts$Taxonomic_distribution ==
                           "R" & results_all_bin_counts$Abundance_distribution %in% c(2, 3), ]
Depth_analysis_2 <-
  results_all_bin_counts[results_all_bin_counts$Taxonomic_distribution ==
                           "R" &
                           results_all_bin_counts$Abundance_distribution == 2, ] # Only using ED 2
Depth_analysis_3 <-
  results_all_bin_counts[results_all_bin_counts$Taxonomic_distribution ==
                           "R" &
                           results_all_bin_counts$Abundance_distribution == 3, ] # Only using ED 3

Depth_results <-
  Depth_analysis %>% group_by(Depth, Abundance_distribution, Pipeline) %>% summarise(
    Average_bins = mean(Number_of_bins),
    Percentage_taxa_found = (mean(Taxonomies_found) / 42 * 100),
    Average_FP = mean(Number_FP)
  ) %>% arrange(Pipeline)
Depth_results_2 <-
  Depth_analysis_2 %>% group_by(Depth, Pipeline) %>% summarise(
    Average_bins = mean(Number_of_bins),
    Percentage_taxa_found = (mean(Taxonomies_found) / 42 * 100),
    Average_FP = mean(Number_FP)
  ) %>% arrange(Pipeline)
Depth_results_3 <-
  Depth_analysis_3 %>% group_by(Depth, Pipeline) %>% summarise(
    Average_bins = mean(Number_of_bins),
    Percentage_taxa_found = (mean(Taxonomies_found) / 42 * 100),
    Average_FP = mean(Number_FP)
  ) %>% arrange(Pipeline)

#write.table(Depth_results,file="Depth_results.tsv",sep = "\t")


##### Evaluate Evenness influence ####

Evenness_analysis <-
  results_all_bin_counts[results_all_bin_counts$Taxonomic_distribution %in% c("C", "N", "V") &
                           results_all_bin_counts$Abundance_distribution %in% c(2, 3, 4, 5) &
                           results_all_bin_counts$Depth == 60, ]

Evenness_results <-
  Evenness_analysis %>% group_by(Taxonomic_distribution, Abundance_distribution, Pipeline) %>% summarise(
    Average_bins = mean(Number_of_bins),
    Percentage_taxa_found = (mean(Taxonomies_found) / 42 * 100),
    Average_FP = mean(Number_FP)
  ) %>% arrange(Pipeline)

#write.table(Evenness_results, file = "Evenness_results.tsv", sep = "\t")


##### Bin quality assessment ####

checkm_table_a <- read.table("01_qual_analysis/CheckM_comm_data/checkm_output_A.tsv",header=T,sep="\t")
checkm_table_b <- read.table("01_qual_analysis/CheckM_comm_data/checkm_output_B.tsv",header=T,sep="\t")
checkm_table_c <- read.table("01_qual_analysis/CheckM_comm_data/checkm_output_C.tsv",header=T,sep="\t")

mock_bins_a <- checkm_table_a[checkm_table_a$Completeness - (checkm_table_a$Contamination*5)>=50,]
mock_bins_b <- checkm_table_b[checkm_table_b$Completeness - (checkm_table_b$Contamination*5)>=50,]
mock_bins_c <- checkm_table_c[checkm_table_c$Completeness - (checkm_table_c$Contamination*5)>=50,]

quality_bins <- rbind(mock_bins_a,mock_bins_b,mock_bins_c)

write.table(quality_bins,file="01_qual_analysis/01_Quality_bins_checkm_R.tsv",sep="\t")


#### Create new tables for each pipeline ####

Pipeline_8k <- quality_bins[quality_bins$Pipeline=="8K",]
Pipeline_DT <- quality_bins[quality_bins$Pipeline=="DT",]
Pipeline_MM <- quality_bins[quality_bins$Pipeline=="MM",]

#### Plots for completeness ####

par(mfrow=c(3,1))

list_pipelines <-
  list("8K" = Pipeline_8k, "DT" = Pipeline_DT, "MM" = Pipeline_MM)

#By Taxonomic relatedness

for (i in 1:length(list_pipelines)) {
  pdf(
    paste0(
      "01_qual_analysis/02_boxplot_comp_contam/boxplot_comp_contam_taxonomy_",
      names(list_pipelines[i]),
      ".pdf"
    )
  )
  boxplot(
    list_pipelines[[i]]$Completeness ~ list_pipelines[[i]]$Taxonomic_distribution,
    ylab = "Completeness",
    xlab = "Taxonomy_profile",
    main = names(list_pipelines[i])
  )
  boxplot(
    list_pipelines[[i]]$Contamination ~ list_pipelines[[i]]$Taxonomic_distribution,
    ylab = "Contamination",
    xlab = "Taxonomy_profile",
    main = names(list_pipelines[i])
  )
  dev.off()
  
}

# Run 1-by-1
boxplot(Pipeline_8k$Completeness ~ Pipeline_8k$Taxonomic_distribution)
boxplot(Pipeline_DT$Completeness ~ Pipeline_DT$Taxonomic_distribution)
boxplot(Pipeline_MM$Completeness ~ Pipeline_MM$Taxonomic_distribution)


##### By Sequencing depth ####

for (i in 1:length(list_pipelines)) {
  pdf(
    paste0(
      "01_qual_analysis/02_boxplot_comp_contam/boxplot_comp_contam_depth_",
      names(list_pipelines[i]),
      ".pdf"
    )
  )
  boxplot(
    list_pipelines[[i]]$Completeness ~ list_pipelines[[i]]$Depth,
    ylab = "Completeness",
    xlab = "Sequencing depth",
    main = names(list_pipelines[i])
  )
  boxplot(
    list_pipelines[[i]]$Contamination ~ list_pipelines[[i]]$Depth,
    ylab = "Contamination",
    xlab = "Sequencing depth",
    main = names(list_pipelines[i])
  )
  dev.off()
  
}

boxplot(Pipeline_8k$Completeness ~ Pipeline_8k$Depth)
boxplot(Pipeline_DT$Completeness ~ Pipeline_DT$Depth)
boxplot(Pipeline_MM$Completeness ~ Pipeline_MM$Depth)

##### By Evenness (abundance) distribution ####

for (i in 1:length(list_pipelines)) {
  pdf(
    paste0(
      "01_qual_analysis/02_boxplot_comp_contam/boxplot_comp_contam_evenness_",
      names(list_pipelines[i]),
      ".pdf"
    )
  )
  boxplot(
    list_pipelines[[i]]$Completeness ~ list_pipelines[[i]]$Abundance_distribution,
    ylab = "Completeness",
    xlab = "Evenness profile",
    main = names(list_pipelines[i])
  )
  boxplot(
    list_pipelines[[i]]$Contamination ~ list_pipelines[[i]]$Abundance_distribution,
    ylab = "Contamination",
    xlab = "Evenness profile",
    main = names(list_pipelines[i])
  )
  dev.off()
  
}

boxplot(Pipeline_8k$Completeness ~ Pipeline_8k$Abundance_distribution)
boxplot(Pipeline_DT$Completeness ~ Pipeline_DT$Abundance_distribution)
boxplot(Pipeline_MM$Completeness ~ Pipeline_MM$Abundance_distribution)
