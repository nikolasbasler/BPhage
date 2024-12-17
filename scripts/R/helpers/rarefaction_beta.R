#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Beta ####
# Full set
beta_start <- Sys.time()
print("Starting beta 1 of 4 with full set.")

met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_dist <- list()
beta_plot_list <- list()
for (tlvl in taxlevels) {
  beta_dist[[tlvl]] <- ordination(df = phage_ab[[tlvl]], 
                                  meta_vars = met_v, 
                                  min_seq = min_seq_count,
                                  df_lengths = phage_lengths[[tlvl]])
  
  beta_plot_list[[tlvl]] <- beta_plot(ordination_list = beta_dist[[tlvl]]$ord_list, 
                                      meta_vars = met_v,
                                      mapped_reads = count_stats$ratios)
}

# Core or not
print("Starting beta 2 of 4 core or not.")

met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_dist_core_or_not <- list()
beta_plot_list_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  for (tlvl in taxlevels) {
    beta_dist_core_or_not[[core_or_not]][[tlvl]] <- ordination(df = phage_ab_core_or_not[[core_or_not]][[tlvl]], 
                                                               meta_vars = met_v, 
                                                               min_seq = min_seq_count,
                                                               df_lengths = phage_lengths[[tlvl]])
    
    beta_plot_list_core_or_not[[core_or_not]][[tlvl]] <- beta_plot(ordination_list = beta_dist_core_or_not[[core_or_not]][[tlvl]]$ord_list, 
                                                                   meta_vars = met_v,
                                                                   mapped_reads = count_stats_core_or_not[[core_or_not]]$ratios)
  }
}

# Absolute counts
print("Starting beta 3 of 4 with absolute counts, full set.")

met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_abs_dist <- list()
beta_abs_plot_list <- list()
for (tlvl in taxlevels) {
  beta_abs_dist[[tlvl]] <- ordination(df = phage_load[[tlvl]],
                                      meta_vars = met_v,
                                      absolute_values = TRUE)
  beta_abs_plot_list[[tlvl]] <- beta_plot(beta_abs_dist[[tlvl]]$ord_list,
                                          meta_vars = met_v,
                                          mapped_reads = count_stats$ratios)
  
}

# Absolute counts core or not
print("Starting beta 4 of 4 with absolute counts core or not.")

met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_abs_dist_core_or_not <- list()
beta_abs_plot_list_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  for (tlvl in taxlevels) {
    beta_abs_dist_core_or_not[[core_or_not]][[tlvl]] <- ordination(df = phage_load_core_or_not[[core_or_not]][[tlvl]],
                                                                   meta_vars = met_v,
                                                                   absolute_values = TRUE)
    beta_abs_plot_list_core_or_not[[core_or_not]][[tlvl]] <- beta_plot(beta_abs_dist[[tlvl]]$ord_list,
                                                                       meta_vars = met_v,
                                                                       mapped_reads = count_stats_core_or_not[[core_or_not]]$ratios)
    
  }
}

beta_end <- Sys.time()
beta_end - beta_start