#0.2) Load libraries

library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)

full_mag_counts <- read.csv("03_pairwise_analysis/00_MAG_counts_full.csv", header = TRUE, sep = ",")

#### Perform Pairwise analyses (Evenness) ####

# Function to perform analysis
analyze_pipeline_eve <- function(pipeline_condition) {
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3"),
           Taxonomy.distribution == "R",
           Pipeline == pipeline_condition)
  
  # Function to process data and perform tests
  process_and_test <- function(data, value_column, file_suffix) {
    data_selected <- data %>%
      arrange(Evenness.distribution) %>%
      select(Evenness.distribution, Sequencing.Depth, all_of(value_column))
    
    summary_data <- data_selected %>%
      group_by(Sequencing.Depth, Evenness.distribution) %>%
      summarise(
        n = n(),
        mean = mean(get(value_column)),
        sd = sd(get(value_column)),
        .groups = 'drop'
      )
    
    formula <- as.formula(paste(value_column, "~ Sequencing.Depth"))
    stat.test <- data_selected %>%
      group_by(Evenness.distribution) %>%
      rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
      select(-.y., -statistic, -df)
    
    write.table(stat.test, file = paste0("03_pairwise_analysis/eve/", file_suffix, "_", pipeline_condition, ".csv"))
    
    return(list("summary_data" = summary_data, "stat.test" = stat.test))
  }
  
  # Function to create plots
  create_plot <- function(data, value_column, plot_file, stat.test) {
    #pdf(plot_file)
    myplot <- ggboxplot(
      data,
      x = "Sequencing.Depth",
      y = value_column,
      fill = "Sequencing.Depth",
      palette = "npg",
      legend = "none"
    ) +
      facet_wrap(~Evenness.distribution)
    
    stat.test <- stat.test$stat.test %>% rstatix::add_xy_position(x = "Sequencing.Depth")
    final_plot <- myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p")
    # Save the plot to a PDF file
    ggplot2::ggsave(plot_file, final_plot, device = "pdf")
    #dev.off()
  }
  
  # MAG numbers
  mag_results <- process_and_test(dat, "Number.of.MAGs", "MAG_counts_Depth_statistics_results")
  create_plot(dat, "Number.of.MAGs", paste0("03_pairwise_analysis/eve/MAG_counts_", pipeline_condition, "_depth_vs_evenness.pdf"), mag_results)
  
  # FP numbers
  fp_results <- process_and_test(dat, "Number.of.False.Positive.Taxonomies", "FP_Depth_statistics_results")
  create_plot(dat, "Number.of.False.Positive.Taxonomies", paste0("03_pairwise_analysis/eve/FP_", pipeline_condition, "_eve_vs_depth.pdf"), fp_results)
  
  # TP numbers
  tp_results <- process_and_test(dat, "Original.taxonomies.found", "TP_Depth_statistics_results")
  create_plot(dat, "Original.taxonomies.found", paste0("03_pairwise_analysis/eve/TP_", pipeline_condition, "_eve_vs_depth.pdf"), tp_results)
  
  # FN numbers
  fn_results <- process_and_test(dat, "Number.of.FalseNegative.Taxonomies", "FN_Depth_statistics_results")
  create_plot(dat, "Number.of.FalseNegative.Taxonomies", paste0("03_pairwise_analysis/eve/FN_Sequencing_depth_", pipeline_condition, "_statistics.pdf"), fn_results)
  
  
}

# New function for specific condition analysis
analyze_specific_condition_eve <- function() {
  dat <- full_mag_counts %>%
    filter(Evenness.distribution == "1",
           Taxonomy.distribution == "O",
           Sequencing.Depth == "60")
  
  process_test_plot <- function(value_column, output_prefix) {
    data_selected <- dat %>%
      arrange(Pipeline) %>%
      select(Pipeline, all_of(value_column))
    
    summary_data <- data_selected %>%
      group_by(Pipeline) %>%
      summarise(
        n = n(),
        mean = mean(get(value_column)),
        sd = sd(get(value_column)),
        .groups = 'drop'
      )
    
    formula <- as.formula(paste(value_column, "~ Pipeline"))
    stat.test <- data_selected %>%
      rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
      select(-.y., -statistic, -df)
    
    write.table(stat.test, file = paste0("03_pairwise_analysis/specific_condition/", output_prefix, "_evenness_statistics_results_pipelines.csv"))
    
    plot_file <- paste0("03_pairwise_analysis/specific_condition/", output_prefix, "_pipelines_60_1_O.pdf")
    create_plot(data_selected, value_column, plot_file, list("stat.test" = stat.test))
  }
  
  create_plot <- function(data, value_column, plot_file, stat.test, width = 8, height = 6) {
    myplot <- ggpubr::ggboxplot(
      data,
      x = "Pipeline",
      y = value_column,
      fill = "Pipeline",
      palette = "npg",
      legend = "none"
    )
    stat.test <- stat.test$stat.test %>% rstatix::add_xy_position(x = "Pipeline")
    final_plot <- myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p")
    
    ggplot2::ggsave(plot_file, final_plot, device = "pdf", width = width, height = height)
  }
  
  # Analyze MAGs, TPs, FPs, and FNs
  process_test_plot("Number.of.MAGs", "MAG_counts")
  process_test_plot("Original.taxonomies.found", "TP")
  process_test_plot("Number.of.False.Positive.Taxonomies", "FP")
  process_test_plot("Number.of.FalseNegative.Taxonomies", "FN")
}

# Call the function with the pipeline condition
analyze_pipeline_eve("MM")
analyze_pipeline_eve("DT")
analyze_pipeline_eve("8K")

#Call specific condition for pipeline
analyze_specific_condition_eve()


#### Perform Pairwise analyses (Taxonomy)#######

analyze_data_tax <- function(pipeline_condition, value_column) {
  # File names are automatically generated based on the value_column and pipeline_condition
  result_file_suffix <- paste0("03_pairwise_analysis/tax/", value_column, "_eve_vs_taxonomy_statistics_results_", pipeline_condition, ".csv")
  plot_file_name <- paste0("03_pairwise_analysis/tax/", value_column, "_", pipeline_condition, "_tax_vs_eve.pdf")
  
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3", "4", "5"),
           Taxonomy.distribution %in% c("C", "N", "V"),
           Sequencing.Depth == "60",
           Pipeline == pipeline_condition) %>%
    arrange(Taxonomy.distribution) %>%
    select(Taxonomy.distribution, Evenness.distribution, all_of(value_column))
  
  summary_data <- dat %>%
    group_by(Taxonomy.distribution, Evenness.distribution) %>%
    summarise(
      n = n(),
      mean = mean(dat[[value_column]]),
      sd = sd(dat[[value_column]]),
      .groups = 'drop'
    )
  
  formula <- as.formula(paste(value_column, "~ Evenness.distribution"))
  stat.test <- dat %>%
    group_by(Taxonomy.distribution) %>%
    rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
    select(-.y., -statistic, -df)
  
  write.table(stat.test, file = result_file_suffix)
  
  ggplot2::ggsave(plot_file_name, {
    ggpubr::ggboxplot(
      dat,
      x = "Evenness.distribution",
      y = value_column,
      fill = "Evenness.distribution",
      palette = "npg",
      legend = "none"
    ) +
      facet_wrap(~ Taxonomy.distribution) +
      ggpubr::stat_pvalue_manual(stat.test %>% rstatix::add_xy_position(x = "Evenness.distribution"), label = "p")
  }, device = "pdf", width = 8, height = 6)
}

# Usage
analyze_data_tax("MM", "Number.of.MAGs")
analyze_data_tax("MM", "Original.taxonomies.found")
analyze_data_tax("MM", "Number.of.False.Positive.Taxonomies")
analyze_data_tax("MM", "Number.of.FalseNegative.Taxonomies")

analyze_data_tax("DT", "Number.of.MAGs")
analyze_data_tax("DT", "Original.taxonomies.found")
analyze_data_tax("DT", "Number.of.False.Positive.Taxonomies")
analyze_data_tax("DT", "Number.of.FalseNegative.Taxonomies")

analyze_data_tax("8K", "Number.of.MAGs")
analyze_data_tax("8K", "Original.taxonomies.found")
analyze_data_tax("8K", "Number.of.False.Positive.Taxonomies")
analyze_data_tax("8K", "Number.of.FalseNegative.Taxonomies")
