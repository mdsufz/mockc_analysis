#0.2) Load libraries

library(dplyr)
library(tidyr)
library(ggplot2)

##################################
##### Pair recovery analysis ####
##################################

#### Check if we always recovered the two closely related species ####

#####

tax_map_merge_edit <- read.csv(file = "04_pair_sp_analysis/03_tax_map_all_edit.csv",
                               row.names = NULL, header = T) #Generated in the initial processing 

tax_map_merge_edit$Species <- gsub("_[A-Z]", "", tax_map_merge_edit$Species)

sample_id <- unique(tax_map_merge_edit$sample)
exp_id <- unique(tax_map_merge_edit$experiment)

######
mprfs_df_format <- read.table(file = "04_pair_sp_analysis/02_rbind_mprfs_all_edit.csv", header = T, row.names = NULL,
                              sep = ",") #Generated in the initial processing 

mprfs_df_format <- unique(mprfs_df_format)

tax_ids <- unique(mprfs_df_format$Tax_dis)

# Found C and V pairs manually
mprfs_df_format_unq_spec <-  mprfs_df_format %>%
  filter(Tax_dis %in% c("V", "C")) %>%
  select(Species) %>%
  unique()

#Open manually edited pairs file
mprfs_species_pairs_vc <- read.csv(file = "04_pair_sp_analysis/04_mockc_final_triplets_table.csv")

mprfs_df_format_vc <-  mprfs_df_format %>%
  filter(Tax_dis %in% c("V", "C"))


recov_pairs_check_list <- list()
pair_precision_df <- data.frame(matrix(ncol = 6, nrow = 0))
pairs_check_list <- list()

for (i in 1:length(sample_id)) {
  
  for (j in 1:length(exp_id)) {
    
    #Select only C and V Tax_dis
    
    filt_df <- tax_map_merge_edit %>%
      filter(experiment == exp_id[j]) %>%
      filter(sample == sample_id[i]) %>%
      filter(Tax_dis %in% c("C", "V")) %>%
      select(Community, Curve, Tax_dis, Depth, rel_abd, coverage, Species)
    
    total_number_recov_bins <- nrow(filt_df)
    
    mprfs_filt <- mprfs_df_format_vc %>%
      filter(sample == sample_id[i]) %>%
      filter(Tax_dis %in% c("C", "V")) %>%
      select(Community, Curve, Tax_dis, Depth, Abd, Species) %>%
      arrange(desc(Abd))
    
    total_number_of_original_species <- nrow(mprfs_filt)
    
    
    if (total_number_of_original_species > 0) {
      
      community_select <- unique(filt_df$Community)
      tax_dis_select <- unique(filt_df$Tax_dis)
      
      original_pairs_df <- mprfs_species_pairs_vc %>%
        filter(Community == community_select)
      
      
      if (tax_dis_select == "V") {
        
        pairs_check_V_df <- data_frame()
        
        for (p in 1:nrow(original_pairs_df)) {
          
          pair_tmp <- original_pairs_df[p, c("Reference","Very.Closely.Related")]
          
          found_pair <- all(pair_tmp %in% filt_df$Species)
          
          row_found <- c(pair_tmp, found_pair)
          names(row_found)[3] <- "found_pair"
          pairs_check_V_df <- rbind.data.frame(pairs_check_V_df, row_found)
          
        }
        
        # Found / Total
        accuracy_pair <- length(which(pairs_check_V_df$found_pair)) / length(pairs_check_V_df$found_pair)
        
        filt_df_select <- filt_df %>%
          select(Species, coverage, rel_abd)
        
        test <- merge.data.frame(pairs_check_V_df,
                                 filt_df_select,
                                 by.x = "Reference",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test) -1, length(test))] <- 
          c("coverage.found.ref", "rel_abd.found.ref")
        
        test <- merge.data.frame(test,
                                 filt_df_select,
                                 by.x = "Very.Closely.Related",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test) -1, length(test))] <- 
          c("coverage.found.Related",
            "rel_abd.found.Related")
        
        test <- merge.data.frame(test,
                                 mprfs_filt %>%
                                   select(Species, Abd),
                                 by.x = "Reference",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test))] <- 
          c("rel_abd.original.ref")
        
        test <- merge.data.frame(test,
                                 mprfs_filt %>%
                                   select(Species, Abd),
                                 by.x = "Very.Closely.Related",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test))] <- 
          c("rel_abd.original.Related")
        
        test <- test %>%
          select(Reference, Very.Closely.Related, everything())
        
        test <- unique(test)
        
        pairs_check_list[[sample_id[i]]][[exp_id[j]]] <- test
        
        
      } else if (tax_dis_select == "C") {
        
        pairs_check_C_df <- data_frame()
        
        for (p in 1:nrow(original_pairs_df)) {
          
          pair_tmp <- original_pairs_df[p, c("Reference","Close.Related")]
          
          found_pair <- all(pair_tmp %in% filt_df$Species)
          
          row_found <- c(pair_tmp, found_pair)
          names(row_found)[3] <- "found_pair"
          
          pairs_check_C_df <- rbind.data.frame(pairs_check_C_df, row_found)
          
        }
        # Found / Total
        accuracy_pair <- length(which(pairs_check_C_df$found_pair)) / length(pairs_check_C_df$found_pair)
        
        filt_df_select <- filt_df %>%
          select(Species, coverage, rel_abd)
        
        test <- merge.data.frame(pairs_check_C_df,
                                 filt_df_select,
                                 by.x = "Reference",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test) -1, length(test))] <- 
          c("coverage.found.ref", "rel_abd.found.ref")
        
        test <- merge.data.frame(test,
                                 filt_df_select,
                                 by.x = "Close.Related",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test) -1, length(test))] <- 
          c("coverage.found.Related",
            "rel_abd.found.Related")
        
        test <- merge.data.frame(test,
                                 mprfs_filt %>%
                                   select(Species, Abd),
                                 by.x = "Reference",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test))] <- 
          c("rel_abd.original.ref")
        
        test <- merge.data.frame(test,
                                 mprfs_filt %>%
                                   select(Species, Abd),
                                 by.x = "Close.Related",
                                 by.y = "Species",
                                 all.x = T)
        
        names(test)[c(length(test))] <- 
          c("rel_abd.original.Related")
        
        test <- test %>%
          select(Reference, Close.Related, everything())
        
        test <- unique(test)
        
        pairs_check_list[[sample_id[i]]][[exp_id[j]]] <- test
        
        
      }
      
      row_acc <- c(unlist(strsplit(x = sample_id[i], "[.]")),
                   exp_id[j],
                   accuracy_pair)
      
      pair_precision_df <- rbind.data.frame(pair_precision_df, row_acc)
      
      #Save to file -  Uncomment to save the files
      # write.csv(pairs_check_list[[sample_id[i]]][[exp_id[j]]],
      #           paste0("04_pair_sp_analysis/10_new_pair_results_per_community/",
      #                  sample_id[i],"_",exp_id[j],".csv"), row.names = F)
      
    }
    names(pair_precision_df) <- c("Community", "Evenness_distribution",
                                  "Taxonomic_distribution", "Sequencing_depth",
                                  "Pipeline", "Found.percentage")
    
  }
}

#Save pair summary to file
write.csv(pair_precision_df, "04_pair_sp_analysis/05_pair_recov_new_triplets_summary.csv", row.names = F)
