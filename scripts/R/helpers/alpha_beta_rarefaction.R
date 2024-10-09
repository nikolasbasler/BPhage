#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Alpha ####
# Full set and core or not
print("Starting alpha 1 of 4 with full set.")
alpha_start <- Sys.time()
met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha <- list()
alpha_core_or_not <- list()
for (tlvl in taxlevels) {
  alpha[[tlvl]] <- alpha_stats(df = phage_ab[[tlvl]], 
                               meta_vars = met_v, 
                               min_seq = min_seq_count,
                               df_lengths = phage_lengths[[tlvl]])
  for (core_or_not in unique(classification$Core)) {
    alpha_core_or_not[[core_or_not]][[tlvl]] <- alpha_stats(df = phage_ab_core_or_not[[core_or_not]][[tlvl]], 
                                                            meta_vars = met_v, 
                                                            min_seq = min_seq_count_core_or_not[[core_or_not]],
                                                            df_lengths = phage_lengths[[tlvl]])
  }
}

# By country and core or not by country
print("Starting alpha 2 of 4 by country and core or not.")
met_v <- c("Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_by_country <- list()
alpha_by_country_core_or_not <- list()
for (countr in levels(metadata$Country)) {
  for (tlvl in taxlevels) {
    count_filt_ab <- phage_ab[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    alpha_by_country[[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                      meta_vars = met_v, 
                                                      min_seq = min_seq_count,
                                                      df_lengths = phage_lengths[[tlvl]])
    for (core_or_not in unique(classification$Core)) {
      count_filt_ab <- phage_ab_core_or_not[[core_or_not]][[tlvl]] %>%
        select(all_of(tlvl), starts_with(countr))
      alpha_by_country_core_or_not[[core_or_not]][[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                                                   meta_vars = met_v, 
                                                                                   min_seq = min_seq_count_core_or_not[[core_or_not]],
                                                                                   df_lengths = phage_lengths[[tlvl]])
    }
  }
}

# Absolute counts full (measured) set and core or not
print("Starting alpha 3 of 4 with absolute counts, full set and core or not.")
met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_abs <- list()
alpha_abs_core_or_not <- list()
for (tlvl in taxlevels) {
  alpha_abs[[tlvl]] <- alpha_stats(df = phage_load[[tlvl]], 
                                   absolut_values = TRUE,
                                   meta_vars = met_v)
  for (core_or_not in unique(classification$Core)) {
    alpha_abs_core_or_not[[core_or_not]][[tlvl]] <- alpha_stats(df = phage_load_core_or_not[[core_or_not]][[tlvl]], 
                                                                absolut_values = TRUE,
                                                                meta_vars = met_v)
  }
}

# Absolute counts by country and core or not by country
print("Starting alpha 4 of 4 with absolute counts cor or not country by country.")

met_v <- c("Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_abs_by_country <- list()
alpha_abs_by_country_core_or_not <- list()
for (countr in levels(metadata$Country)) {
  for (tlvl in taxlevels) {
    count_filt_ab <- phage_load[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    alpha_abs_by_country[[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                          absolut_values = TRUE,
                                                          meta_vars = met_v)
    for (core_or_not in unique(classification$Core)) {
      count_filt_ab <- phage_load_core_or_not[[core_or_not]][[tlvl]] %>%
        select(all_of(tlvl), starts_with(countr))
      alpha_abs_by_country_core_or_not[[core_or_not]][[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                                                       absolut_values = TRUE,
                                                                                       meta_vars = met_v)
    }
  }
}
alpha_end <- Sys.time()
alpha_end - alpha_start
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