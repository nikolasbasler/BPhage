report_stats <- function(df, thresholds) {
  outlist <- list()
  num_of_seqs <- df %>%
    as.data.frame() %>%
    column_to_rownames(colnames(df)[1]) %>%
    t() %>% 
    as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    pivot_longer(-Sample_ID) %>%
    group_by(Sample_ID) %>%
    summarize(n_seq = sum(value)) %>%
    mutate(ratio_to_highest = max(n_seq)/n_seq) %>%
    mutate(gut_part = str_split(Sample_ID, "_", simplify = TRUE)[, 4])
  
  outlist$ratios <- num_of_seqs %>%
      select(-gut_part) %>%
      arrange(n_seq)
  
  outlist$plot_hist <- ggplot(num_of_seqs, aes(x=n_seq)) +
    geom_histogram(binwidth = 100000)
  
  outlist$plot_densities <- ggplot(num_of_seqs, aes(x=n_seq, color=gut_part)) +
    geom_density()
  
  outlist$plot_thresholds <- num_of_seqs %>%
    mutate(gut_part = str_split(Sample_ID, "_", simplify = TRUE)[, 4]) %>%
    group_by(gut_part) %>%
    arrange(n_seq) %>%
    mutate(sample_number = 1:n()) %>%
    ggplot(aes(x=sample_number, y=n_seq, color=gut_part)) +
    geom_point() +
    scale_y_continuous(trans="log10") +
    geom_hline(yintercept = thresholds)
  
  read_stats = tibble()
  for (t in c(0, thresholds)) {
    discarded <- num_of_seqs %>%
      filter(n_seq <= t) %>%
      group_by(gut_part) %>%
      summarise(n_sample = n()) %>%
      as.data.frame() %>%
      column_to_rownames("gut_part") %>%
      t() %>%
      as.data.frame()
    
    r_stats <- num_of_seqs %>%
      filter(n_seq>t) %>%
      summarise(threshold=t,
                lowest = min(n_seq), 
                highest = max(n_seq),
                extreme_ratio = max(ratio_to_highest),
                median_ratio = median(ratio_to_highest),
                quantile_90_ratio = quantile(ratio_to_highest, probs=0.90),
                quantile_95_ratio = quantile(ratio_to_highest, probs=0.95),
                samples_left = length(n_seq),
                with_ratio_above_1000 = sum(ratio_to_highest>1000),
                with_ratio_above_2500 = sum(ratio_to_highest>2500),
                discarded_samples = nrow(num_of_seqs)-samples_left
      ) %>%
      mutate(mid = discarded$mid) %>%
      mutate(ile = discarded$ile) %>%
      mutate(rec = discarded$rec)
    read_stats <- bind_rows(read_stats, r_stats)
  }
  outlist$read_stats <- read_stats
  return(outlist)
}

discards <- function(ratios, min_seq) {
  bees_before_filter <- ratios %>%
    mutate(bees = str_sub(Sample_ID, start=1, end=12)) %>%
    select(bees) %>%
    unique()
  bees_after_filter <- ratios %>%
    filter(n_seq >= min_seq) %>%
    mutate(bees = str_sub(Sample_ID, start=1, end=12)) %>%
    select(bees) %>%
    unique()
  lost_bees <- setdiff(bees_before_filter, bees_after_filter)
  
  discarded <- ratios %>%
    filter(n_seq < min_seq) %>%
    select(Sample_ID) %>%
    unique()
  return(list(discarded=discarded, lost_bees=lost_bees))
}

