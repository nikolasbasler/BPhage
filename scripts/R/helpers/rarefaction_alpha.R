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

# Get detailed results of pair-wise comparison tests:
shannon_and_meta <- list()
shannon_and_meta$all <- alpha$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)
shannon_and_meta$noncore <- alpha_core_or_not$no$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)
shannon_and_meta$core <- alpha_core_or_not$yes$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)

shannon_and_meta_abs <- list()
shannon_and_meta_abs$all <- alpha_abs$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)
shannon_and_meta_abs$noncore <- alpha_abs_core_or_not$no$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)
shannon_and_meta_abs$core <- alpha_abs_core_or_not$yes$Family$table %>%
  left_join(., metadata[c("Sample_ID", "Gut_part", "Country", "Season")], by = "Sample_ID") %>%
  select(Sample_ID, Hill_Shannon, Gut_part, Country, Season)

alpha_family_pwc <- list()
alpha_abs_family_pwc <- list()
for (set in names(shannon_and_meta)) {
  metavars <- c("Gut_part", "Country", "Season")
  for (mv in metavars) {
    alpha_family_pwc[[set]][[mv]] <- pwc_shannon_square(
      s_and_h = shannon_and_meta[[set]],
      meta_var = mv
    )
  }
  metavars <- c("Country", "Season")
  for (mv in metavars) {
    alpha_abs_family_pwc[[set]][[mv]] <- pwc_shannon_square(
      s_and_h = shannon_and_meta_abs[[set]],
      meta_var = mv
    )
  }
}

alpha_end <- Sys.time()
alpha_end - alpha_start