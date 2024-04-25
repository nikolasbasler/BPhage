
###### BETA DIVERSITY
beta_rarified <- function(ab_table, sampling_depth, lengths, seed) {
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
  
  df <- rrarefy(ab_filt, sample=sampling_depth) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(tax) %>%
    inner_join(., lengths, by=tax) %>%
    calc_reads_per_kb(., tax) %>%
    column_to_rownames(tax) %>%
    t() %>%
    vegdist(method="bray") %>% 
    as.matrix() %>% 
    as.data.frame()
  return(df)
}

distance_histogram <- function(dist_df, metadata_var, metadata_value) {
  dist_df %>%
    rownames_to_column("Sample_ID") %>%
    pivot_longer(-Sample_ID, values_to = "distance") %>%
    filter(Sample_ID > name) %>%
    ggplot(aes(x=distance)) +
    geom_histogram(binwidth = 0.01) +
    labs(title = paste0(metadata_var, " - ", metadata_value))
}

rared_ordination = function(df, min_seq, meta_vars, df_lengths, seed) {
  tax <- colnames(df)[1]
  df_t = df %>% 
    filter(.data[[tax]]!="unclassified") %>%
    column_to_rownames(tax) %>% 
    t() %>% 
    as.data.frame()
  
  plan(multisession, workers=4)
  beta_df_list = future_map(1:iterations,
                             ~beta_rarified(ab_table = df_t,
                                             sampling_depth = min_seq,
                                             lengths = df_lengths,
                                             seed = .x)) 
  beta_average_df <- Reduce(`+`, beta_df_list) / length(beta_df_list) 
  beta_average_dist <- as.dist(beta_average_df)

  ord_list <- list()
  ord_list$all$all <- pcoa(beta_average_dist)
  
  distances_plots <- list()
  distances_plots$all$all <- beta_average_dist %>%
    as.matrix() %>% 
    as.data.frame() %>% 
    distance_histogram(metadata_var = "all", metadata_value = "all")
  
  for (m_var in meta_vars) {
    for (m_value in levels(metadata[[m_var]])) {
      filter_vector <- metadata %>%
        filter(.[[m_var]]==m_value) %>%
        select(Sample_ID) %>%
        unlist() %>%
        as.character()
      
        beta_average_filt <- beta_average_dist %>%
        as.matrix() %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample_ID") %>%
        select(Sample_ID, any_of(filter_vector)) %>%
        filter(Sample_ID %in% filter_vector) %>%
        column_to_rownames("Sample_ID")
      
        ord_list[[m_var]][[m_value]] <- pcoa(beta_average_filt)
        
        distances_plots[[m_var]][[m_value]] <- beta_average_filt %>%
          distance_histogram(metadata_var = m_var, metadata_value = m_value)
    }
  }
  return(list(ord_list = ord_list, avg_dist = beta_average_df, dist_hist_list = distances_plots))
}

beta_plot = function(ordination_list, meta_vars, mapped_reads) {
  final_plotlist <- list()
  for (n in names(ordination_list)) {
    for (m in names(ordination_list[[n]])) {
      plot_df <- metadata %>%
        select(Sample_ID, all_of(meta_vars)) %>%
        inner_join(., rownames_to_column(as.data.frame(ordination_list[[n]][[m]]$vectors),"Sample_ID"),
                   by="Sample_ID") %>%
        left_join(., mapped_reads, by = "Sample_ID")
      
      plotlist_for_patch <- list()
      for (meta_v in meta_vars) {
        if (n==meta_v) {
          next
        }
        plotlist_for_patch[[meta_v]] <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=.data[[meta_v]])) +
          labs(x=NULL, y=NULL) +
          geom_point()
      }
      patch <- wrap_plots(plotlist_for_patch)
      
      # This is just to give the different panels one axis label.
      x_axis_lab <- paste0("PCo1 (", round(ordination_list[[n]][[m]]$values$Relative_eig[1]*100,2),"%)")
      y_axis_lab <- paste0("PCo2 (", round(ordination_list[[n]][[m]]$values$Relative_eig[2]*100,2),"%)")
      empty_y_plot <- ggplot(data.frame(l = y_axis_lab, x = 1, y = 1)) +
        geom_text(aes(x, y, label = l), angle = 90) +
        theme_void() +
        coord_cartesian(clip = "off")
      empty_x_plot <- ggplot(data.frame(l = x_axis_lab, x = 1, y = 1)) +
        geom_text(aes(x, y, label = l)) +
        theme_void() +
        coord_cartesian(clip = "off")
      
      final_plotlist[[n]][[m]] <- (empty_y_plot + patch + plot_layout(widths = c(1, 45) )) / empty_x_plot + plot_layout(heights = c(45, 1)) +
        plot_annotation(title=paste0(n,": ",m))
    }
    
    
    if (n=="all") {
      plot_df <- metadata %>%
        select(Sample_ID, all_of(meta_vars)) %>%
        inner_join(., rownames_to_column(as.data.frame(ordination_list$all$all$vectors),"Sample_ID"),
                   by="Sample_ID") %>%
        left_join(., mapped_reads, by = "Sample_ID")
      final_plotlist$all$control <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=n_seq)) +
        geom_point() +
        labs(x=x_axis_lab, y=y_axis_lab) +
        plot_annotation(title="all: control")
    }
  }
  return(final_plotlist)
}
