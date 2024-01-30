# Library
library(fmsb)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(gridGraphics)
library(png)

all_comm_results <-
  read.table("02_All_comms_results_R_means.txt",
             header = T,
             sep = "\t") #File created on initial_data step

all_comm_results$Depth <- as.factor(all_comm_results$Depth)

all_comm_results$Abundance_distribution <-
  as.factor(all_comm_results$Abundance_distribution)

plot_data_depth_abdV <- function(data, numerical_column, abundance_value, plot_title, filename) {
  
  # Group and summarize the data
  data_summary <- data %>%
    filter(Taxonomic_distribution == "R") %>%
    group_by(Depth, Abundance_distribution, Pipeline) %>%
    summarise(mean_value = mean({{numerical_column}}), .groups = 'drop') %>%
    filter(Abundance_distribution == abundance_value) %>%
    ungroup()
  
  data_summary$Depth <- paste0(data_summary$Depth, " mio")
  
  # Pivot the data
  data_pivot <- data_summary %>%
    select(-Abundance_distribution) %>%
    pivot_wider(values_from = mean_value, names_from = Depth) %>%
    as.data.frame()
  
  row.names(data_pivot) <- data_pivot$Pipeline
  data_pivot <- data_pivot %>% select(-Pipeline)
  
  # Calculate max value for the plot and adjust axis labels
  max_value <- max(data_summary$mean_value, na.rm = TRUE)
  max_plot_value <- ceiling(max_value / 10) * 10 # Round up to the nearest 10
  axis_labels <- seq(0, max_plot_value, by = max_plot_value / 4)
  
  # Add min and max lines
  data_pivot <- rbind(rep(max_plot_value, ncol(data_pivot)), rep(0, ncol(data_pivot)), data_pivot)
  
  # Define colors (customize as needed)
  colors_border <- c(rgb(230/255,74/255,51/255,0.9),
                     rgb(76/255,186/255,212/255,0.9),
                     rgb(0/255,163/255,137/255,0.9))
  colors_in <- c(rgb(230/255,74/255,51/255,0.4),
                 rgb(76/255,186/255,212/255,0.4),
                 rgb(0/255,163/255,137/255,0.4))
  
  # Plot
  png(filename, width = 800, height = 600)
  radarchart(data_pivot, axistype=1,
             pcol=colors_border, pfcol=colors_in, plwd=4, plty=1,
             cex.main = 2,
             cglcol="grey", cglty=1, axislabcol="black",calcex=2, caxislabels=axis_labels, cglwd=0.8,
             vlcex=2.5, title = plot_title)
  # Prepare the legend as a separate plot
legend(x=1, y=1.3, legend = rownames(data_pivot[-c(1,2),]), bty = "n", pch=20, col=colors_border, text.col = "black", cex=3, pt.cex=5)


  dev.off()
  
  }


# Usage

plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Taxonomies_found, 
                     abundance_value = "2", "True Positive on Logarithmic decay",
                     filename = "tp_log_spider.png")

plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Number_FP, 
                     abundance_value = "2", "Number of False Positives on Logarithmic decay",
                     filename = "fp_log_spider.png")

plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Number_of_bins, 
                     abundance_value = "2", "Number of MAGs on Logarithmic decay",
                     filename = "nmags_log_spider.png")


img1 <- readPNG("nmags_log_spider.png")
img2 <- readPNG("tp_log_spider.png")
img3 <- readPNG("fp_log_spider.png")


plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Taxonomies_found, 
                     abundance_value = "3", "True Positive on Exponential decay",
                     filename = "tp_exp_spider.png")

plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Number_FP, 
                     abundance_value = "3", "Number of False Positives on Exponential decay",
                     filename = "fp_exp_spider.png")

plot_data_depth_abdV(data = all_comm_results, 
                     numerical_column = Number_of_bins, 
                     abundance_value = "3", "Number of MAGs on Exponential decay",
                     filename = "nmags_exp_spider.png")

img4 <- readPNG("nmags_exp_spider.png")
img5 <- readPNG("tp_exp_spider.png")
img6 <- readPNG("fp_exp_spider.png")


grid.arrange(rasterGrob(img1), rasterGrob(img2), rasterGrob(img3),
             rasterGrob(img4), rasterGrob(img5), rasterGrob(img6), nrow = 2)

