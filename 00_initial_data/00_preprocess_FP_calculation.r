#0.1) Set working dir

#0.2) Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)

#0.3) Load taxa data from bins

all_bin_tax <- read.table(file = "00_initial_data/01_all_generated_bins_taxonomy_edit.csv", sep = "\t",
                          header = T) #Results from GTDBtk
names(all_bin_tax)[1] <- "bin_id"


#0.4) Load mprfs outputs systematically (Mock community files)

files <- dir(recursive=TRUE,
             full.names=TRUE,
             pattern="\\.mprf$")


#Start dataframe with one example

mprfs_list <- list()

sample_name <- sub(".*mprfs/", "", files[1])
sample_name <- sub(".mprf*", "", sample_name)
sample_name <- gsub(pattern = "_", replacement = ".", sample_name)

mprfs_list[[sample_name]] <- read.table(file = files[1])
names(mprfs_list[[sample_name]]) <- c("Abd", "name", "Species")
mprfs_list[[sample_name]] <- mprfs_list[[sample_name]][,c(1,3)]

mprfs_list[[sample_name]]$sample <- sample_name

mprfs_df <- mprfs_list[[sample_name]]

for (i in 2:length(files)) {
  
  sample_name <- sub(".*mprfs/", "", files[i])
  sample_name <- sub(".mprf*", "", sample_name)
  sample_name <- gsub(pattern = "_", replacement = ".", sample_name)
  
  mprfs_list[[sample_name]] <- read.table(file = files[i])
  names(mprfs_list[[sample_name]]) <- c("Abd", "name", "Species")
  mprfs_list[[sample_name]] <- mprfs_list[[sample_name]][,c(1,3)]
  
  mprfs_list[[sample_name]]$sample <- sample_name
  
  mprfs_df <- rbind.data.frame(mprfs_df, mprfs_list[[sample_name]])
  
}

#Save to file

write.csv(mprfs_df, file = "00_initial_data/02_rbind_mprfs_all.csv", row.names = F)

#Manually edit the mockc_data/02_rbind_mprfs_all.csv Species name

mprfs_df_format <- read.table(file = "00_initial_data/02_rbind_mprfs_all_edit.csv", header = T, row.names = NULL,
                              sep = ",")

mprfs_df_format <- unique(mprfs_df_format)


# 0.5) Load mapped table

#Coverage

map_reduced_cov_list <- read.table(file = "00_initial_data/map_tables/map_reduced_coverage_list.tsv",
                                   sep = " ")

names(map_reduced_cov_list) <- c("bin_id", "sample", "coverage")


#relative abundance
map_reduced_relabd_list <- read.table(file = "00_initial_data/map_tables/map_reduced_relative_abundance_list.tsv",
                                      sep = " ")

names(map_reduced_relabd_list) <- c("bin_id", "sample", "rel_abd")


#### 1) Merge gtdb bin taxa with cov list #####

relabd_cov_merge <- merge.data.frame(x = map_reduced_cov_list,
                                     y = map_reduced_relabd_list,
                                     by = "bin_id")


tax_map_merge <- merge.data.frame(x = relabd_cov_merge,
                                  y = all_bin_tax,
                                  by = "bin_id")


#### 1.1) Save gtdb merged files to manually extract 8k, DT.. ####

write.csv(tax_map_merge, file = "00_initial_data/03_tax_map_all.csv", row.names = F)

#Read data
tax_map_merge_edit <- read.csv(file = "00_initial_data/03_tax_map_all_edit.csv",
                               row.names = NULL, header = T)

tax_map_merge_edit$Species <- gsub("_[A-Z]", "", tax_map_merge_edit$Species)


##### 2) Generate TP / FP table ####

sample_id <- unique(tax_map_merge_edit$sample)
exp_id <- unique(tax_map_merge_edit$experiment)

found_species_df <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(found_species_df) <- c("sample",
                                "experiment",
                                "total_recov_bins",
                                "total_orig_bins",
                                "precision",
                                "recall",
                                "f1_score",
                                "n_TP",
                                "n_FP",
                                "avg_TP_cov",
                                "avg_FP_cov",
                                "min_tp_cov",
                                "max_fp_cov")

merge_tax_mprfs_rank_list <- list()
mprf_list <- list()
filt_df_list <- list()

for (i in 1:length(sample_id)) {
  for (j in 1:length(exp_id)) {
    
    filt_df <- tax_map_merge_edit %>%
      filter(experiment == exp_id[j]) %>%
      filter(sample == sample_id[i]) %>%
      select(rel_abd, coverage, Species)
    
    total_number_recov_bins <- nrow(filt_df)
    
    mprfs_filt <- mprfs_df_format %>%
      filter(sample == sample_id[i]) %>%
      select(Abd, Species) %>%
      arrange(desc(Abd))
    
    total_number_of_original_species <- nrow(mprfs_filt)
    
    if (nrow(filt_df)>0) {
      
      filt_df <- filt_df %>% 
        aggregate(rel_abd ~ Species + coverage,. , sum) %>%
        arrange(desc(rel_abd))
      
      filt_df <- filt_df %>% arrange(desc(coverage))
      
      merge_tax_mprfs_filt <- merge.data.frame(x = filt_df,
                                               y = mprfs_filt,
                                               by = "Species") %>% arrange(desc(rel_abd))
      
      
    }else{
      
      merge_tax_mprfs_filt <- merge.data.frame(x = filt_df,
                                               y = mprfs_filt,
                                               by = "Species") %>% arrange(desc(rel_abd))
      
    }
    
    mprf_list[[sample_id[i]]] <- mprfs_filt
    
    #Calculate metrics
    #TRUE POSITIVE
    # Found in the recov and present in the original
    filt_df$true_positive <- filt_df$Species %in% mprfs_filt$Species
    true_positive <- sum(filt_df$Species %in% mprfs_filt$Species)
    
    avg_tp_cov <- filt_df %>%
      filter(true_positive == TRUE) %>%
      summarise(mean(coverage)) %>%
      as.numeric(.) %>%
      replace_na(0)
    
    min_tp_cov <- filt_df %>%
      filter(true_positive == TRUE) %>%
      summarise(min(coverage)) %>%
      mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
      as.numeric(.) %>%
      replace_na(0)
    
    min_tp_specie <- filt_df %>%
      filter(coverage == min_tp_cov) %>%
      select(Species)
    
    #TRUE NEGATIVE
    # not found in the recov, not present in the original
    
    #FALSE POSITIVE
    # Found in the recov, but not present in the original
    filt_df$false_positive <- !(filt_df$Species %in% mprfs_filt$Species)
    false_positive <- sum(!(filt_df$Species %in% mprfs_filt$Species))
    
    avg_fp_cov <- filt_df %>%
      filter(false_positive == TRUE) %>%
      summarise(mean(coverage)) %>%
      as.numeric(.) %>%
      replace_na(0)
    
    max_fp_cov <- filt_df %>%
      filter(false_positive == TRUE) %>%
      summarise(max(coverage)) %>%
      mutate_all(function(x) ifelse(is.infinite(x), 0, x)) %>%
      as.numeric(.) %>%
      replace_na(0)
    
    
    #FALSE NEGATIVE
    # not found in the recov, but present in the original
    false_negative <- sum(!(mprfs_filt$Species %in% filt_df$Species))
    
    #calculate metrics
    precision <- true_positive/(true_positive+false_positive)
    
    recall <- true_positive/(true_positive+false_negative)
    
    f1_score <- 2*precision*recall/(precision+recall)
    
    #
    community <- tax_map_merge_edit %>%
      filter(sample == sample_id[i]) %>%
      select(Community) %>% unique(.) %>%
      as.character(.)
    
    if (nrow(filt_df) > 0) {
      
      filt_df$Community <- community
      filt_df$sample <- sample_id[i]
      filt_df$pipeline <- exp_id[j]
      
      #Save FP Taxa
      # write.csv(filt_df %>%
      #       filter(false_positive == TRUE) %>%
      #         select(Species, sample, pipeline, Community, coverage, rel_abd),
      #     paste0("11_false_positives_df/FP_",
      #            sample_id[i],"_",exp_id[j],".csv"), row.names = F)
      
    }
    
    filt_df_list[[sample_id[i]]][[exp_id[j]]] <- filt_df
    
    
    #Check if they are the found species keep same abd order
    
    mprfs_filt_onlyfound <- mprfs_filt %>%
      filter(Species %in% filt_df$Species)
    
    select_row <- data.frame(sample_id[i],
                             exp_id[j],
                             total_number_recov_bins,
                             total_number_of_original_species,
                             precision,
                             recall,
                             f1_score,
                             true_positive,
                             false_positive,
                             avg_tp_cov,
                             avg_fp_cov,
                             min_tp_cov,
                             max_fp_cov)
    
    #print(select_row)
    found_species_df <- rbind(found_species_df, select_row)
    
  }
}

#Boxplot of AVG cov FP and AVG cov TP

cov_tp_fp_df <- rbind.data.frame(found_species_df %>%
                                   select(exp_id.j., avg_tp_cov) %>%
                                   rename(., COV = avg_tp_cov) %>%
                                   mutate(class = "TP"),
                                 found_species_df %>%
                                   select(exp_id.j., avg_fp_cov) %>%
                                   rename(., COV = avg_fp_cov) %>%
                                   mutate(class = "FP"))

ggplot(cov_tp_fp_df, aes(x=exp_id.j., y=COV, fill=class)) + 
  geom_boxplot()


found_species_df_split <- tidyr::separate(found_species_df, col = sample_id.i.,
                                          into = c("community", "abd_distribution", "tax_distribution", "seq_depth"), sep = "\\.")

write.csv(found_species_df_split, file = "00_initial_data/04_All_comms_results_R.csv")

#After some small manual edit, the found_species_df_split becomes All_comms_results_R.csv (Supp Table 2)

