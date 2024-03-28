
# calc_tpm <- function(abtable, length_df) {
#   abtable %>%
#     rownames_to_column("contig") %>%
#     inner_join(., length_df, by="contig") %>% 
#     mutate(across(-contig, ~./length_kb)) %>% 
#     select(-length_kb) %>% 
#     mutate(across(-contig, ~./sum(.)))
# }
# calc_reads_per_kb <- function(abtable, length_df) {
#   abtable %>%
#     rownames_to_column("contig") %>%
#     inner_join(., length_df, by="contig") %>%
#     mutate(across(-contig, ~./length_kb)) %>% 
#     select(-length_kb) %>%
#     mutate(across(-contig, ~round(.)))
# }

calc_tpm <- function(abtable, level, lengths_df) {
  abtable %>%
    inner_join(., lengths_df, by=level) %>%
    mutate(across(-all_of(level), ~./length_kb)) %>% 
    select(-length_kb) %>% 
    mutate(across(-all_of(level), ~./sum(.)))
}


calc_reads_per_kb <- function(abtable_with_lengths, level) {
  abtable_with_lengths %>%
    mutate(across(-all_of(level), ~./length_kb)) %>% 
    select(-length_kb) # %>%
    # mutate(across(-all_of(level), ~round(.)))
}


tax_lvl_per_kb = function(tax_level, contig_lengths, ab_table, classif ) {
  sample_names=colnames(ab_table)
  ab_table <- ab_table %>%
    rownames_to_column("contig") %>%
    inner_join(., contig_lengths, by="contig") 

  classif <- classif  %>%
    rownames_to_column("contig") %>%
    mutate_all(~ ifelse(. == "", "unclassified", .))
  
  inner_join(ab_table, classif, by="contig") %>% 
    group_by(.data[[tax_level]]) %>%
    summarize_at(vars(all_of(c(sample_names,"length_kb"))), sum) %>%
    mutate(across(-all_of(tax_level), ~./length_kb)) %>% 
    select(-length_kb) %>%
    mutate(across(-all_of(tax_level), ~round(.))) %>%
    select_if(~ !(is.numeric(.) && all(. == 0)))
}

tax_sum = function(tax_level, ab_table, classif) {
  sample_names <- colnames(ab_table)[2:length(ab_table)]
  classif <- classif  %>%
    mutate(contig = as.character(contig)) %>%
    mutate_all(~ ifelse(. == "", "unclassified", .))
  
  inner_join(ab_table, classif, by="contig") %>%
    group_by(.data[[tax_level]]) %>%
    summarize_at(vars(all_of(sample_names)), sum)
}

tax_lengths <- function(tax_level, classif) {
  classif %>%
    mutate(contig = as.character(contig)) %>%
    group_by(.data[[tax_level]]) %>%
    summarise(length_kb = sum(length_kb)) %>%
    mutate_all(~ ifelse(. == "", "unclassified", .))
}

hostg_filter <- function(hg, ab_table, classif) {
  classif %>%
    filter(Host_group==hg) %>%
    inner_join(., ab_table, by="contig") %>%
    select(colnames(ab_table)) %>% 
    select_if(~ any(. > 0))
}

