alpha_rarefied = function(ab_table, sampling_depth, lengths, seed) {
  set.seed(seed)
  ab_filt <- ab_table %>%
    filter(rowSums(.) >= sampling_depth)
  if (nrow(ab_table) > nrow(ab_filt)) {
    warning(nrow(ab_table)-nrow(ab_filt), 
            " samples discarded because they had fewer than ", sampling_depth, 
            " observations: ", 
            paste(setdiff(rownames(ab_table), rownames(ab_filt)), collapse=", ")
            )
  }
  tax <- colnames(lengths)[1]

  df = rrarefy(ab_filt, sample=sampling_depth) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(tax) %>%
    inner_join(., lengths, by=tax) %>%
    calc_reads_per_kb(., tax) %>%
    column_to_rownames(tax) %>%
    t() %>% 
    as_tibble(rownames="sample") %>%
    group_by(sample) %>%
    pivot_longer(-sample) %>%
    summarize(Richness = specnumber(value),
              # evenness = shannon_idx/log(observed_species), # This is Pilou's evenness
              shannon_idx = diversity(value, index="shannon"),
              Hill_Shannon = exp(shannon_idx),
              Hill_Simpson = diversity(value, index="invsimpson")
    ) %>%
    select(-shannon_idx) %>%
    as.data.frame() %>%
    column_to_rownames("sample")
  df
}

alpha_stats = function(df, meta_vars, min_seq = NA, df_lengths = NA, absolut_values = FALSE) {
  meta_vars <- factor(meta_vars, levels = c("Gut_part", "Country", "Season", "Health"))
  
  tax <- colnames(df)[1]
  df_t = df %>% 
    filter(.data[[tax]]!="Unclassified") %>%
    column_to_rownames(tax) %>% 
    t() %>% 
    as.data.frame()
  
  if (absolut_values) {
    alpha_tbl <- df_t %>%
      as_tibble(rownames="Sample_ID") %>%
      pivot_longer(-Sample_ID) %>%
      group_by(Sample_ID) %>%
      summarize(Richness = specnumber(value),
                # evenness = shannon_idx/log(observed_species), # This is Pilou's evenness
                shannon_idx = diversity(value, index="shannon"),
                Hill_Shannon = exp(shannon_idx),
                Hill_Simpson = diversity(value, index="invsimpson")
      ) %>%
      select(-shannon_idx)
    
  } else {
    plan(multisession, workers=4)
    alpha_df_list = future_map(1:iterations,
                               ~alpha_rarefied(ab_table = df_t,
                                               sampling_depth = min_seq,
                                               lengths = df_lengths,
                                               seed = .x)) 
    alpha_average_df = Reduce(`+`, alpha_df_list) / length(alpha_df_list)
    alpha_tbl = as_tibble(alpha_average_df, rownames = "Sample_ID")
  }

  pre_plot_tibble <- metadata %>%
    select("Sample_ID", all_of(meta_vars)) %>%
    inner_join(., alpha_tbl, by="Sample_ID") %>%
    pivot_longer(all_of(meta_vars), names_to = "meta_variable", values_to = "meta_value") %>%
    pivot_longer(c(Richness, Hill_Shannon, Hill_Simpson), names_to = "metric") %>%
    mutate(metric = factor(metric, levels=c("Richness", "Hill_Shannon", "Hill_Simpson")))
  
  box_fill_colors <- list(Gut_part = "#ef8f01", 
                          Country = "#8B4513", 
                          Season = "#FFA07A", 
                          Health = "white")

  plot_tibble <- pre_plot_tibble %>%
    group_by(metric, meta_variable) %>%
    pairwise_wilcox_test(
      value ~ meta_value,
      p.adjust.method = "BH"
    ) %>%
    group_by(metric, meta_variable) %>%
    summarise(
      cld = list({
        pvec <- setNames(p.adj, paste(group1, group2, sep = "-"))
        multcompLetters(pvec, threshold = 0.05)$Letters
      }), .groups = "drop"
    ) %>%
    unnest_longer(cld, values_to = "letter", indices_to = "meta_value") %>%
    full_join(pre_plot_tibble, ., by = c("meta_variable", "meta_value", "metric")) %>%
    mutate(meta_value = factor(meta_value, levels = levels(pre_plot_tibble$meta_value)))

  kruskal_results <- plot_tibble %>%
    group_by(meta_variable, metric) %>%
    summarize(pvalue = kruskal.test(value~meta_value)$p.value,
              test_stat = kruskal.test(value~meta_value)$statistic,
              deg_freedom = kruskal.test(value~meta_value)$parameter,
              .groups="drop")

  # This is necessary, because geom_pwc() doesn't work well with faceted plots...
  panels <- list()
  number_of_metrics <- length(levels(plot_tibble$metric))
  number_of_meta_vars <- length(meta_vars)
  total_number_of_plots <- number_of_metrics * number_of_meta_vars
  for (metr in levels(plot_tibble$metric)) {
    # for (meta_v in levels(meta_vars)) {
    for (meta_v in meta_vars) {
      krusk <- kruskal_results %>%
        filter(metric==metr &  meta_variable==meta_v)
      stats_text = paste0("KW test ",
                          # round(krusk$test_stat, digits=2), " .",
                          "p-value: ",
                          round(krusk$pvalue, digits=2))
      
      annotation_df <- plot_tibble %>%
        filter(metric == metr,
               meta_variable == meta_v) %>%
        group_by(meta_value) %>%
        mutate(sample_count = n()) %>%
        ungroup() %>%
        mutate(letter_y = max(value)*1.12) %>%
        select(-c(Sample_ID, value)) %>%
        distinct()
      
      p <- plot_tibble %>%
        filter(metric == metr) %>%
        filter(meta_variable == meta_v) %>%
        ggplot(aes(x = meta_value, y = value, fill = meta_variable)) +
        geom_boxplot() +
        labs(x = NULL, y = NULL) +
        scale_fill_manual(values = box_fill_colors[[meta_v]]) +
        guides(fill="none") +
        
        # Puts sample number at the bottom of the boxes. Makes the plot very busy...
        # geom_text(data = annotation_df, 
        #           mapping = aes(x = meta_value, y = 0, label = sample_count),
        #           inherit.aes = FALSE,
        #           size = 3,
        #           vjust = 1
        # ) +
        
        theme_minimal()
      
      if (krusk$pvalue > 0.05) {
        annotation_df <- annotation_df %>%
          mutate(letter = "")
      }

      p <- p + 
        geom_label(data = annotation_df, 
                  mapping = aes(x = meta_value, y = letter_y, label = letter),
                  inherit.aes = FALSE,
                  fill        = "white",
                  label.size  = 0,
                  label.padding = unit(0.35, "lines")
                  # size = 5
                  )
      if ((length(panels)+1) %% number_of_meta_vars == 1) {
        p <- p + labs(y=metr)
      }
      if (total_number_of_plots - length(panels) <= number_of_meta_vars) {
        p <- p + labs(x=meta_v)
      }
      panels[[length(panels)+1]] <- p
    }
  }
  patch <- wrap_plots(panels, ncol=length(meta_vars), axes = "collect_y")
  return(list(plot=patch, table=alpha_tbl, kruskal = kruskal_results, single_plots = panels))
}

pwc_shannon_square <- function(s_and_h, meta_var) {
  
  groups <- levels(factor(s_and_h[[meta_var]]))  
  pwc <- rstatix::pairwise_wilcox_test(s_and_h, formula = as.formula(paste("Hill_Shannon ~", meta_var)), p.adjust.method = "BH") %>%
    select(group1, group2, statistic, p.adj)
  
  square <- pwc %>%
    select(-p.adj) %>%
    complete(group1 = groups, group2 = groups, fill = list(statistic = NA)) %>%
    mutate(
      group1 = factor(group1, levels = groups),
      group2 = factor(group2, levels = groups)
    ) %>%
    arrange(group1, group2) %>%
    left_join(., pwc[c("group1", "group2", "p.adj")], by = join_by(group1 == group2, group2 == group1)) %>%
    mutate(value = case_when(
      is.finite(statistic) & is.na(p.adj) ~ statistic,
      is.finite(p.adj) & is.na(statistic) ~ p.adj,
      .default = NA
    )) %>%
    select(group1, group2, value) %>%
    pivot_wider(names_from = group2) %>%
    rename(group = group1)
  return(square)
}
