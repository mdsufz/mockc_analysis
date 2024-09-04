#0.2) Load libraries

library(tidyr);packageVersion("tidyr")
library(ggplot2);packageVersion("ggplot2")
library(dplyr);packageVersion("dplyr")
library(ggpubr);packageVersion("ggpubr")
library(rstatix);packageVersion("rstatix")
library(car);packageVersion("car")  # for Levene's test
library(emmeans);packageVersion("emmeans")  # for Tukey HSD

full_mag_counts <- read.csv("03_pairwise_analysis/00_MAG_counts_full.csv", header = TRUE, sep = ",")

#### Perform Pairwise analyses (Evenness) ####

# Function to perform analysis
analyze_pipeline_eve <- function(pipeline_condition) {
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3"),
           Taxonomy.distribution == "R",
           Pipeline == pipeline_condition)
  
  process_and_test <- function(data, value_column, file_suffix) {
    data_selected <- data %>%
      arrange(Evenness.distribution) %>%
      select(Evenness.distribution, Sequencing.Depth, all_of(value_column)) %>%
      mutate(
        Sequencing.Depth = as.factor(Sequencing.Depth),
        Evenness.distribution = as.factor(Evenness.distribution)
      )
    
    summary_data <- data_selected %>%
      group_by(Sequencing.Depth, Evenness.distribution) %>%
      summarise(
        n = n(),
        mean = mean(get(value_column)),
        sd = sd(get(value_column)),
        .groups = 'drop'
      )
    
    # T-test with Bonferroni correction
    formula <- as.formula(paste(value_column, "~ Sequencing.Depth"))
    stat.test <- data_selected %>%
      group_by(Evenness.distribution) %>%
      rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
      select(-.y., -statistic, -df)
    
    # ANOVA
    anova_formula <- as.formula(paste(value_column, "~ Sequencing.Depth * Evenness.distribution"))
    anova_result <- aov(anova_formula, data = data_selected)
    anova_summary <- summary(anova_result)
    
    # Tukey HSD
    tukey_result <- TukeyHSD(anova_result)
    
    # Save results
    write.table(stat.test, file = paste0("03_pairwise_analysis/eve/ttest_", file_suffix, "_", pipeline_condition, ".csv"), row.names = F, sep = ",")
    write.table(anova_summary[[1]], file = paste0("03_pairwise_analysis/eve/anova_", file_suffix, "_", pipeline_condition, ".csv"), sep = ",")
    write.table(tukey_result$`Sequencing.Depth:Evenness.distribution`, file = paste0("03_pairwise_analysis/eve/tukey_", file_suffix, "_", pipeline_condition, ".csv"), sep = ",")
    
    return(list("summary_data" = summary_data, "stat.test" = stat.test, 
                "anova_summary" = anova_summary, "tukey_result" = tukey_result))
  }
  
  create_plot <- function(data, value_column, plot_file, results) {
    myplot <- ggboxplot(
      data,
      x = "Sequencing.Depth",
      y = value_column,
      fill = "Sequencing.Depth",
      palette = "npg",
      legend = "none"
    ) +
      facet_wrap(~Evenness.distribution)
    
    stat.test <- results$stat.test %>% rstatix::add_xy_position(x = "Sequencing.Depth")
    final_plot <- myplot + 
      ggpubr::stat_pvalue_manual(stat.test, label = "p")
    
    print(final_plot)
    
    ggplot2::ggsave(plot_file, final_plot, device = "pdf")
  }
  
  # Process each metric
  metrics <- c("Number.of.MAGs", "Number.of.False.Positive.Taxonomies", "Original.taxonomies.found", "Number.of.FalseNegative.Taxonomies")
  file_suffixes <- c("MAG_counts_Depth_statistics_results", "FP_Depth_statistics_results", "TP_Depth_statistics_results", "FN_Depth_statistics_results")
  
  for (i in 1:length(metrics)) {
    results <- process_and_test(dat, metrics[i], file_suffixes[i])
    create_plot(dat, metrics[i], paste0("03_pairwise_analysis/eve/", metrics[i], "_", pipeline_condition, "_depth_vs_evenness.pdf"), results)
  }
}

# New function for specific condition analysis
analyze_specific_condition_eve <- function() {
  dat <- full_mag_counts %>%
    filter(Evenness.distribution == "1",
           Taxonomy.distribution == "O",
           Sequencing.Depth == "60") %>%
    mutate(Pipeline = as.factor(Pipeline))
  
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
    
    # T-test with Bonferroni correction
    formula <- as.formula(paste(value_column, "~ Pipeline"))
    stat.test <- data_selected %>%
      rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
      select(-.y., -statistic, -df)
    
    # One-way ANOVA
    anova_result <- aov(formula, data = data_selected)
    anova_summary <- summary(anova_result)
    
    # Tukey HSD
    tukey_result <- TukeyHSD(anova_result)
    
    # Save results
    write.table(summary_data, file = paste0("03_pairwise_analysis/specific_condition/summary_", output_prefix, "_evenness_statistics_results_pipelines.csv"), row.names = F, sep = ",")
    write.table(stat.test, file = paste0("03_pairwise_analysis/specific_condition/ttest_", output_prefix, "_evenness_statistics_results_pipelines.csv"), row.names = F, sep = ",")
    write.table(anova_summary[[1]], file = paste0("03_pairwise_analysis/specific_condition/anova_", output_prefix, "_evenness_statistics_results_pipelines.csv"), sep = ",")
    write.table(tukey_result$Pipeline, file = paste0("03_pairwise_analysis/specific_condition/tukey_", output_prefix, "_evenness_statistics_results_pipelines.csv"), sep = ",")
    
    # Create plot
    plot_file <- paste0("03_pairwise_analysis/specific_condition/", output_prefix, "_pipelines_60_1_O.pdf")
    plot <- create_plot(data_selected, value_column, list("stat.test" = stat.test, "anova_summary" = anova_summary, "tukey_result" = tukey_result))
    
    # Save plot
    ggplot2::ggsave(plot_file, plot, device = "pdf", width = 10, height = 8)
    
    return(list(summary_data = summary_data, 
                t_test = stat.test, 
                anova = anova_summary, 
                tukey = tukey_result))
  }
  
  create_plot <- function(data, value_column, results) {
    myplot <- ggpubr::ggboxplot(
      data,
      x = "Pipeline",
      y = value_column,
      fill = "Pipeline",
      palette = "npg",
      legend = "none"
    )
    
    # Add t-test p-values
    stat.test <- results$stat.test %>% rstatix::add_xy_position(x = "Pipeline")
    myplot <- myplot + ggpubr::stat_pvalue_manual(stat.test, label = "p.adj")
    
    # Add ANOVA and Tukey information
    anova_p <- results$anova_summary[[1]]["Pipeline", "Pr(>F)"]
    tukey_sig <- sum(results$tukey_result$Pipeline[,"p adj"] < 0.05)
    

    return(myplot)
  }
  
  # Analyze MAGs, TPs, FPs, and FNs
  results <- list(
    MAG_counts = process_test_plot("Number.of.MAGs", "MAG_counts"),
    TP = process_test_plot("Original.taxonomies.found", "TP"),
    FP = process_test_plot("Number.of.False.Positive.Taxonomies", "FP"),
    FN = process_test_plot("Number.of.FalseNegative.Taxonomies", "FN")
  )
  
  return(results)
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
  result_file_suffix <- paste0("03_pairwise_analysis/tax/", value_column, "_eve_vs_taxonomy_statistics_results_", pipeline_condition)
  plot_file_name <- paste0("03_pairwise_analysis/tax/", value_column, "_", pipeline_condition, "_tax_vs_eve.pdf")
  
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3", "4", "5"),
           Taxonomy.distribution %in% c("C", "N", "V"),
           Sequencing.Depth == "60",
           Pipeline == pipeline_condition) %>%
    arrange(Taxonomy.distribution) %>%
    select(Taxonomy.distribution, Evenness.distribution, all_of(value_column)) %>%
    mutate(
      Taxonomy.distribution = as.factor(Taxonomy.distribution),
      Evenness.distribution = as.factor(Evenness.distribution)
    )
  
  summary_data <- dat %>%
    group_by(Taxonomy.distribution, Evenness.distribution) %>%
    summarise(
      n = n(),
      mean = mean(get(value_column)),
      sd = sd(get(value_column)),
      .groups = 'drop'
    )
  
  # T-test with Bonferroni correction
  formula <- as.formula(paste(value_column, "~ Evenness.distribution"))
  stat.test <- dat %>%
    group_by(Taxonomy.distribution) %>%
    rstatix::t_test(formula, p.adjust.method = "bonferroni") %>%
    select(-.y., -statistic, -df)
  
  # Two-way ANOVA
  anova_formula <- as.formula(paste(value_column, "~ Taxonomy.distribution * Evenness.distribution"))
  anova_result <- aov(anova_formula, data = dat)
  anova_summary <- summary(anova_result)
  
  # Tukey HSD for Evenness distribution within each Taxonomy distribution
  tukey_results <- dat %>%
    group_by(Taxonomy.distribution) %>%
    do(tukey = TukeyHSD(aov(as.formula(paste(value_column, "~ Evenness.distribution")), data = .)))
  
  # Save results
  write.table(stat.test, file = paste0(result_file_suffix, "_ttest.csv"), row.names = F, sep = ",")
  write.table(anova_summary[[1]], file = paste0(result_file_suffix, "_anova.csv"), sep = ",")
  
  # Save Tukey results for each Taxonomy distribution
  for (tax in levels(dat$Taxonomy.distribution)) {
    write.table(tukey_results$tukey[[which(tukey_results$Taxonomy.distribution == tax)]]$Evenness.distribution, 
                file = paste0(result_file_suffix, "_tukey_", tax, ".csv"), sep = ",")
  }
  
  # Create plot
  plot <- ggpubr::ggboxplot(
    dat,
    x = "Evenness.distribution",
    y = value_column,
    fill = "Evenness.distribution",
    palette = "npg",
    legend = "none"
  ) +
    facet_wrap(~ Taxonomy.distribution) +
    ggpubr::stat_pvalue_manual(stat.test %>% rstatix::add_xy_position(x = "Evenness.distribution"), label = "p")
  
  ggplot2::ggsave(plot_file_name, plot, device = "pdf", width = 10, height = 8)
  
  return(list(summary_data = summary_data, 
              t_test = stat.test, 
              anova = anova_summary, 
              tukey = tukey_results))
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


############## Depth 60MM ###############
analyze_data_pipeline_depth_taxo <- function(value_column, pipeline) {
  result_file_suffix <- paste0("03_pairwise_analysis/tax_eve/", value_column, "_", pipeline, "_comparison_results")
  plot_file_name <- paste0("03_pairwise_analysis/tax_eve/", value_column, "_", pipeline, "_comparison_results.pdf")
  
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3", "4", "5"),
           Taxonomy.distribution %in% c("C", "V", "N"),
           Sequencing.Depth %in% c("60"),
           Pipeline == pipeline) %>%
    arrange(Sequencing.Depth, Evenness.distribution) %>%
    select(Evenness.distribution, Taxonomy.distribution, all_of(value_column)) %>%
    mutate(
      Evenness.distribution = as.factor(Evenness.distribution),
      Taxonomy.distribution = as.factor(Taxonomy.distribution)
    )
  
  summary_data <- dat %>%
    group_by(Taxonomy.distribution, Evenness.distribution) %>%
    summarise(
      n = n(),
      mean = mean(get(value_column)),
      sd = sd(get(value_column)),
      .groups = 'drop'
    )
  
  # Pairwise t-tests with Bonferroni correction
  stat.test <- dat %>%
    group_by(Evenness.distribution) %>%
    rstatix::t_test(as.formula(paste(value_column, "~ Taxonomy.distribution")), 
                    p.adjust.method = "bonferroni") %>%
    select(-`.y.`, -statistic, -df)
  
  # Two-way ANOVA
  anova_formula <- as.formula(paste(value_column, "~ Taxonomy.distribution * Evenness.distribution"))
  anova_result <- aov(anova_formula, data = dat)
  anova_summary <- summary(anova_result)
  
  # Tukey HSD
  tukey_result <- TukeyHSD(anova_result)
  
  # Save results
  write.csv(summary_data, file = paste0(result_file_suffix, "_summary.csv"), row.names = FALSE)
  write.csv(stat.test, file = paste0(result_file_suffix, "_ttest.csv"), row.names = FALSE)
  write.csv(anova_summary[[1]], file = paste0(result_file_suffix, "_anova.csv"))
  write.csv(tukey_result$Taxonomy.distribution, file = paste0(result_file_suffix, "_tukey_taxonomy.csv"))
  write.csv(tukey_result$Evenness.distribution, file = paste0(result_file_suffix, "_tukey_evenness.csv"))
  write.csv(tukey_result$`Taxonomy.distribution:Evenness.distribution`, file = paste0(result_file_suffix, "_tukey_interaction.csv"))
  
  # Create plot
  plot <- ggboxplot(
    dat,
    x = "Taxonomy.distribution",
    y = value_column,
    fill = "Taxonomy.distribution",
    palette = "npg",
    facet.by = c("Evenness.distribution"),
    legend = "none"
  ) +
    stat_pvalue_manual(
      stat.test %>% add_xy_position(x = "Taxonomy.distribution"),
      label = "p.adj",
      step.increase = 0.1
    ) +
    labs(
      title = paste("Comparison of", value_column, "across Taxonomy.distribution"),
      x = "Taxonomy.distribution",
      y = value_column
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(plot_file_name, plot, device = "pdf", width = 15, height = 10)
  
  return(list(summary_data = summary_data, 
              t_test = stat.test,
              anova = anova_summary,
              tukey = tukey_result))
}

# Usage
analyze_data_pipeline_depth_taxo("Number.of.MAGs", "8K")
analyze_data_pipeline_depth_taxo("Original.taxonomies.found", "8K")
analyze_data_pipeline_depth_taxo("Number.of.False.Positive.Taxonomies", "8K")
analyze_data_pipeline_depth_taxo("Number.of.FalseNegative.Taxonomies", "8K")

analyze_data_pipeline_depth_taxo("Number.of.MAGs", "DT")
analyze_data_pipeline_depth_taxo("Original.taxonomies.found", "DT")
analyze_data_pipeline_depth_taxo("Number.of.False.Positive.Taxonomies", "DT")
analyze_data_pipeline_depth_taxo("Number.of.FalseNegative.Taxonomies", "DT")

analyze_data_pipeline_depth_taxo("Number.of.MAGs", "MM")
analyze_data_pipeline_depth_taxo("Original.taxonomies.found", "MM")
analyze_data_pipeline_depth_taxo("Number.of.False.Positive.Taxonomies", "MM")
analyze_data_pipeline_depth_taxo("Number.of.FalseNegative.Taxonomies", "MM")

############## Pipelines ###############
#Depth
analyze_data_pipeline_itself_depth <- function(value_column) {
  result_file_suffix <- paste0("03_pairwise_analysis/pipelines/", value_column, "_comparison_results")
  plot_file_name <- paste0("03_pairwise_analysis/pipelines/", value_column, "_comparison_results.pdf")
  
  dat <- full_mag_counts %>%
    filter(Evenness.distribution %in% c("2", "3"),
           Taxonomy.distribution %in% c("R"),
           #Sequencing.Depth %in% c("60"),
    ) %>%
    arrange(Sequencing.Depth, Evenness.distribution) %>%
    select(Evenness.distribution, Sequencing.Depth, Pipeline, all_of(value_column)) %>%
    mutate(
      Evenness.distribution = as.factor(Evenness.distribution),
      Sequencing.Depth = as.factor(Sequencing.Depth)
    )
  
  summary_data <- dat %>%
    group_by(Sequencing.Depth, Evenness.distribution) %>%
    summarise(
      n = n(),
      mean = mean(get(value_column)),
      sd = sd(get(value_column)),
      .groups = 'drop'
    )
  
  # Pairwise t-tests with Bonferroni correction
  stat.test <- dat %>%
    group_by(Evenness.distribution, Sequencing.Depth) %>%
    rstatix::t_test(as.formula(paste(value_column, "~ Pipeline")), 
                    p.adjust.method = "bonferroni") %>%
    select(-`.y.`, -statistic, -df)
  
  # One ANOVA
  anova_formula <- as.formula(paste(value_column, "~ Pipeline * Evenness.distribution * Sequencing.Depth"))
  anova_result <- aov(anova_formula, data = dat)
  anova_summary <- summary(anova_result)
  
  # Tukey HSD
  tukey_result <- TukeyHSD(anova_result)
  
  # Save results
  write.csv(summary_data, file = paste0(result_file_suffix, "_summary.csv"), row.names = FALSE)
  write.csv(stat.test, file = paste0(result_file_suffix, "_ttest.csv"), row.names = FALSE)
  write.csv(anova_summary[[1]], file = paste0(result_file_suffix, "_anova.csv"))
  write.csv(tukey_result$`Pipeline:Evenness.distribution:Sequencing.Depth`, file = paste0(result_file_suffix, "_tukey_interaction.csv"))
  
  return(list(summary_data = summary_data, 
              t_test = stat.test,
              anova = anova_summary,
              tukey = tukey_result))
}

# Usage
analyze_data_pipeline_itself_depth("Number.of.MAGs")
analyze_data_pipeline_itself_depth("Original.taxonomies.found")
analyze_data_pipeline_itself_depth("Number.of.False.Positive.Taxonomies")
analyze_data_pipeline_itself_depth("Number.of.FalseNegative.Taxonomies")

#Taxa
analyze_data_pipeline_itself_taxa <- function(value_column) {
  result_file_suffix <- paste0("03_pairwise_analysis/pipelines/", value_column, "_taxa_comparison_results")
  plot_file_name <- paste0("03_pairwise_analysis/pipelines/", value_column, "_taxa_comparison_results.pdf")
  
  dat <- full_mag_counts %>%
    filter(#Evenness.distribution %in% c("2", "3"),
      Taxonomy.distribution %in% c("C", "N", "V"),
      Sequencing.Depth %in% c("60"),
    ) %>%
    arrange(Taxonomy.distribution, Evenness.distribution) %>%
    select(Evenness.distribution, Taxonomy.distribution, Pipeline, all_of(value_column)) %>%
    mutate(
      Evenness.distribution = as.factor(Evenness.distribution),
      Taxonomy.distribution = as.factor(Taxonomy.distribution)
    )
  
  summary_data <- dat %>%
    group_by(Taxonomy.distribution, Evenness.distribution) %>%
    summarise(
      n = n(),
      mean = mean(get(value_column)),
      sd = sd(get(value_column)),
      .groups = 'drop'
    )
  
  # Pairwise t-tests with Bonferroni correction
  stat.test <- dat %>%
    group_by(Evenness.distribution, Taxonomy.distribution) %>%
    rstatix::t_test(as.formula(paste(value_column, "~ Pipeline")), 
                    p.adjust.method = "bonferroni") %>%
    select(-`.y.`, -statistic, -df)
  
  # One ANOVA
  anova_formula <- as.formula(paste(value_column, "~ Pipeline * Evenness.distribution * Taxonomy.distribution"))
  anova_result <- aov(anova_formula, data = dat)
  anova_summary <- summary(anova_result)
  
  # Tukey HSD
  tukey_result <- TukeyHSD(anova_result)
  
  # Save results
  write.csv(summary_data, file = paste0(result_file_suffix, "_summary.csv"), row.names = FALSE)
  write.csv(stat.test, file = paste0(result_file_suffix, "_ttest.csv"), row.names = FALSE)
  write.csv(anova_summary[[1]], file = paste0(result_file_suffix, "_anova.csv"))
  write.csv(tukey_result$`Pipeline:Evenness.distribution:Taxonomy.distribution`, file = paste0(result_file_suffix, "_tukey_interaction.csv"))
  
  return(list(summary_data = summary_data, 
              t_test = stat.test,
              anova = anova_summary,
              tukey = tukey_result))
}

# Usage
analyze_data_pipeline_itself_taxa("Number.of.MAGs")
analyze_data_pipeline_itself_taxa("Original.taxonomies.found")
analyze_data_pipeline_itself_taxa("Number.of.False.Positive.Taxonomies")
analyze_data_pipeline_itself_taxa("Number.of.FalseNegative.Taxonomies")
