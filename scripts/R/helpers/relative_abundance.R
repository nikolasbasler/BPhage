
host_colors <- c("Bifidobacterium" = "#FFDAB9", "Lactobacillus" = "#FFA07A",
                 "Snodgrassella" = "#FFC300", "Bombilactobacillus" = "#ef8f01",
                 "Gilliamella" = "#D2691E", "Frischella" = "#8B4513",
                 "Bartonella" = "#1C3A3A", "Bombella" = "black",
                 "other" = "#555555", "unknown" = "lightgrey")
# prev_colors <- c("#FFEB99", "#EBD688", "#D9BF77", "#C6A866", "#B39155", "#997A3C", "#806600", "black")
prev_colors <- c("#FFDAB9", "#FFA07A", "#FFC300", "#ef8f01", "#D2691E", "#8B4513", "#1C3A3A", "black")


# rel_abund_taxlvl = function(df, tl, meta_vars) {
#   threshold_for_other = 0.005
#   # Todo: Set threshold dynamically by meta_var and tl
#   
#   plotlist <- list()
#   for (m_var in meta_vars) {
#     metadata_filt <- metadata %>%
#       select(Sample_ID, all_of(m_var))
#     
#     plotlist[[m_var]] <- df %>% 
#       column_to_rownames("name") %>% 
#       t() %>%
#       as.data.frame() %>%
#       rownames_to_column("Sample_ID") %>%
#       inner_join(., metadata_filt, by="Sample_ID") %>% 
#       select(-Sample_ID) %>%
#       pivot_longer(-all_of(m_var)) %>%
#       group_by(.data[[m_var]]) %>%
#       mutate(value = value/sum(value)) %>% # Here I'm adding normalized counts (TPM) from different samples. Feels strange, but if I would go back to the read counts instead, then samples with more seq depth could swamp more shallow samples.
#       mutate(name = 
#                ifelse(value < threshold_for_other, 
#                       paste0("other (<", round(threshold_for_other*100, digits = 1),"%)"), 
#                       name)) %>%
#       ungroup() %>%
#       group_by(.data[[m_var]], name) %>%
#       mutate(added_value = sum(value)) %>%  # This part is
#       ungroup() %>%                         # only to merge
#       select(-value) %>%                    # identically-
#       distinct() %>%                        # colored shapes.
#       ggplot(aes(x=.data[[m_var]], y=added_value, fill=name)) +
#       geom_col() +
#       ylab("Cumulated TPM") +
#       labs(fill = tl)
#   }
#   return(plotlist)
# }

###### HEATMAP
contig_heatmap <- function(df, classif) {
  df %>% 
    pivot_longer(-contig) %>%
    inner_join(., metadata, by=join_by(name==Sample_ID)) %>%
    inner_join(., classif, by="contig") %>%
    arrange(Phylum, Class, Order, Family, Genus, Species) %>%
    mutate(contig_and_tax = paste0(contig, " [", Family, ";", Genus, ";",Species, "]")) %>%
    mutate(contig_and_tax = factor(contig_and_tax, levels=unique(contig_and_tax))) %>%
    mutate(value = ifelse(value==0, NA, value)) %>%
    rename(sample = name,
           TPM=value) %>%
    mutate(sample = factor(sample, levels = levels(metadata$Sample_ID))) %>%
    ggplot(aes(x=sample, y=contig_and_tax, fill=TPM)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_row(~Country+Season, scales="free_x", space="free") +
    scale_fill_gradient(limits = c(0, 1))
}

make_heatmap = function(df, tl, meta_vars) {
  threshold_for_removal = 0.005
  plotlist <- list()
  for (m_var in meta_vars) {
    # Use bray on rows and columns
    # cite: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-45#MOESM3
    # https://rdrr.io/cran/NeatMap/man/heatmap1.html

    metadata_filt <- metadata %>%
      select(Sample_ID, all_of(m_var))
    
    if (m_var == "Sample_ID") {
      condensed <- df %>% 
        column_to_rownames("name") %>% 
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Sample_ID") %>%
        pivot_longer(-Sample_ID) %>% 
        rename(variable = Sample_ID) 
    } else {
      # condensed <- df %>%
      #   column_to_rownames("name") %>% 
      #   t() %>%
      #   as.data.frame() %>%
      #   rownames_to_column("Sample_ID") %>%
      #   inner_join(., metadata_filt, by="Sample_ID") %>%
      #   select(-Sample_ID) %>%
      #   pivot_longer(-all_of(m_var)) %>%
      #   rename(variable = all_of(m_var)) %>%
      #   group_by(variable) %>%
      #   mutate(value = sum(value)) %>% # Here I'm adding normalized counts (TPM) from different samples. Feels strange, but if I would go back to the read counts instead, then samples with more seq depth could swamp more shallow samples.
      #   group_by(name) %>%
      #   filter(any(value > threshold_for_removal)) %>%
      #   ungroup() %>%
      #   group_by(variable, name) %>%
      #   mutate(value = sum(value)) %>%        # This part is only
      #   ungroup() %>%                         # to merge identically-
      #   distinct() #%>%                        # colored shapes.
      #   # group_by(variable) %>%
      #   # mutate(value = value/sum(value))
      
      condensed <- df %>%
        column_to_rownames("name") %>% 
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Sample_ID") %>%
        # inner_join(., metadata_filt, by="Sample_ID") %>% select(Sample_ID, `Deformed wing virus`, Health) %>% group_by(Health) %>% mutate(DWV=sum(`Deformed wing virus`)) %>% View() 
        inner_join(., metadata_filt, by="Sample_ID") %>%
        select(-Sample_ID) %>%
        pivot_longer(-all_of(m_var)) %>%
        rename(variable = all_of(m_var)) %>%
        group_by(name, variable) %>% 
        mutate(value = sum(value)/n()) %>% # Here I'm adding normalized counts (TPM) from different samples. Feels strange, but if I would go back to the read counts instead, then samples with more seq depth could swamp more shallow samples.
        ungroup() %>%
        distinct() %>% 
        group_by(name) %>%
        filter(any(value >= threshold_for_removal))

    }
    condensed_ord <- condensed %>%
      pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
      column_to_rownames("variable") %>%
      vegdist() %>%
      pcoa()
    # samples_angles = apply(condensed_ord$vectors[,1:2],1,function(x){atan2(x[1],x[2])})
    # samples_order = names(samples_angles)[order(samples_angles)] # Circular order
    samples_order = names(condensed_ord$vectors[,1])[order(condensed_ord$vectors[,1])] # Order based on 1st component
    # condensed_ord_t <- condensed %>%
    #   pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
    #   column_to_rownames("variable") %>%
    #   t() %>%
    #   vegdist() %>%
    #   pcoa()
    # taxa_angles = apply(condensed_ord_t$vectors[,1:2],1,function(x){atan2(x[1],x[2])})
    # taxa_order = names(taxa_angles)[order(taxa_angles)] # Circular order
    # taxa_order = names(condensed_ord_t$vectors[,1])[order(condensed_ord_t$vectors[,1])] # Order based on 1st component
    taxa_order = df$name[order(df$name, decreasing = TRUE)] # Alphabetic taxon order
    plotlist[[m_var]] <- condensed %>%
      # mutate(value= ifelse(value==0, 0, sqrt(value))) %>% # Square root transformation.
      # mutate(value= ifelse(value==0, 0, log(value))) %>% # Log transformation.
      mutate(name = factor(name, levels=taxa_order),
             variable = factor(variable, levels=samples_order)) %>%
      ggplot(aes(x=variable, y=name, fill=value)) +
      geom_tile() +
      labs(x=NULL, y=NULL, fill="avg. TPM") +
      ggtitle(paste0(tl, "/",m_var," (>", round(threshold_for_removal*100, digits = 1),"%)"))
    
    
    # Base R heatmap with hierarchical clustering for comparison
    # df_log = df %>%
    #   column_to_rownames("name") %>%
    #   log() %>%
    #   mutate_all(~ifelse(. == -Inf, 0, .))
    # df_log %>%
    #   as.matrix() %>%
    #   heatmap()
  }
  return(plotlist)
}

#### PHYLOSEQ HEAT MAP

# data("GlobalPatterns")
# sample_data(GlobalPatterns)

phylo_heat_map <- function(ab_table, id_table) {
  otu_mat <- ab_table %>% 
    column_to_rownames("name") %>%
    as.matrix()
  tax_mat <- id_table %>% 
  column_to_rownames("name") %>%
    as.matrix()
  
  OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  SAM <- sample_data(metadata[colnames(otu_mat),])
  
  ps <- phyloseq(OTU, TAX, SAM)
  
  sub_ps <- subset_taxa(ps, Phylum!="")
  
  plot_heatmap(sub_ps, sample.label=Country, taxa.label="Phylum")
  
  
  # plot_bar(ps, fill = "Family")
}

average_tpm_bar_plot <- function(tpm_table, tl, hg, meta_vars, title_prefix="", threshold_for_other=0.01, hg_or_core = "") {
  label_for_other <- paste0("other (<", round(threshold_for_other*100, digits = 1),"%)")
  tax <- colnames(tpm_table)[1]
  plot_list <- list()
  tible_list <- list()
  for (m_var in meta_vars) {
    metadata_filt <- metadata %>%
      select(Sample_ID, all_of(m_var))
    
    phages_in_group <- phage_abundance %>%
      pivot_longer(-contig, names_to = "Sample_ID") %>%
      mutate(present = ifelse(value > 0, TRUE, FALSE)) %>%
      filter(present) %>%
      left_join(., metadata_filt, by="Sample_ID") %>%
      select(all_of(c("contig", m_var))) %>%
      distinct() %>%
      group_by(.data[[m_var]]) %>%
      summarise(in_group = n())
    
    tible_list[[m_var]] <- tpm_table %>%
      rename(group = all_of(tax)) %>%
      pivot_longer(-group, names_to = "Sample_ID") %>%
      left_join(., metadata_filt, by="Sample_ID") %>%
      mutate(Sample_ID = factor(Sample_ID, levels=sample_order)) %>% 
      filter(!is.na(value)) %>%
      group_by(.data[[m_var]], group) %>%
      summarise(mean_tpm = mean(value), .groups="drop") %>%
      mutate(group = 
               ifelse(mean_tpm < threshold_for_other, 
                      label_for_other, 
                      group)) %>%
      group_by(.data[[m_var]], group) %>%
      mutate(mean_tpm = sum(mean_tpm)) %>%  # This part is only to
      ungroup() %>%                         # merge identically-
      distinct() %>%                        # colored shapes.
      left_join(., phages_in_group, by = m_var)
    
    if (m_var == "Sample_ID") {
      tible_list[[m_var]] <- tible_list[[m_var]] %>%
        left_join(., metadata, by = "Sample_ID") %>%
        select(Sample_ID, group, mean_tpm, in_group, Country, Hive_ID, Season, Gut_part)
    }
    if (tax == c("Prevalence_Countries")) {
      tible_list[[m_var]] <- tible_list[[m_var]] %>%
        mutate(group = as.factor(group))
    }
    if (tl == "Host_group") {
      tible_list[[m_var]] <- tible_list[[m_var]] %>%
        mutate(group = factor(group, levels = rev(c("Bifidobacterium", "Lactobacillus", "Snodgrassella",
                                                   "Bombilactobacillus", "Gilliamella", "Frischella",
                                                   "Bartonella", "Bombella", "other", "unknown"))))
    }
    
    if (tl == "Prevalence") {
      max_val <- max(as.integer(tible_list[[m_var]]$group))
      ticks <- 1:max_val %>% 
        quantile(probs = c(0.33, 0.66)) %>% 
        round()
      
      # # BURN AFTER READING!
      # order_by_phages_in_group <- tible_list[[m_var]] %>%
      #   arrange(in_group) %>%
      #   select(Sample_ID) %>%
      #   distinct()
      # tible_list[[m_var]] <- tible_list[[m_var]] %>%
      #   mutate(Sample_ID = factor(Sample_ID, levels = order_by_phages_in_group$Sample_ID))
      
      if (m_var == "Sample_ID") {
        country_plots <- list()
        for (country in unique(tible_list[[m_var]]$Country)) {
          country_plots[[country]] <- tible_list[[m_var]] %>%
            filter(Country == country) %>% 
            arrange(desc(group)) %>%
            ggplot(aes(x=.data[[m_var]], y=mean_tpm, fill = group)) +
            geom_col() +
            ggtitle(paste0(country, title_prefix, hg_or_core, ": \"", hg,"\"")) +
            labs(fill=tl) +
            scale_fill_gradient(low = "#F0F0F0", high = "black",
                                breaks = c(1, ticks, max_val),
                                labels = c(1, ticks, max_val)) +
            guides(fill = guide_colourbar(reverse = TRUE)) +
            labs(fill=tl) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank()) +
            theme(axis.text.x = element_text(angle = 45, hjust=1)) +
            scale_x_discrete(expand = c(0.025, 0)) +
            geom_text(aes(label=in_group, y = 1.01, vjust = 0, angle = 45)) +
            facet_wrap(Season~Hive_ID, scales = "free_x", nrow = 1)
          if (tax == "Prevalence_Countries") {
            country_plots[[country]] <- tible_list[[m_var]] %>%
              filter(Country == country) %>% 
              ggplot(aes(x=.data[[m_var]], y=mean_tpm, fill = as.factor(group))) +
              geom_col() +
              ggtitle(paste0(country, title_prefix, hg_or_core, ": \"", hg,"\"")) +
              labs(fill=tl) +
              scale_fill_manual(values = prev_colors) +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank()) +
              theme(axis.text.x = element_text(angle = 45, hjust=1)) +
              scale_x_discrete(expand = c(0.025, 0)) +
              geom_text(aes(label=in_group, y = 1.01, vjust = 0, angle = 45)) +
              facet_wrap(Season~Hive_ID, scales = "free_x", nrow = 1)
          }
        }
        plot_list[[m_var]] <- wrap_plots(country_plots, ncol = 1)
        
      } else {
        plot_list[[m_var]] <- tible_list[[m_var]] %>%
          arrange(desc(group)) %>%
          ggplot(aes(x=.data[[m_var]], y=mean_tpm, fill = group)) +
          geom_col() +
          ggtitle(paste0(title_prefix, hg_or_core, ": \"", hg,"\"")) +
          labs(fill=tl) +
          theme_minimal() +
          theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
      }
      if (tax == c("Prevalence_Countries")) {
          plot_list[[m_var]] <- plot_list[[m_var]] +
            scale_fill_manual(values = prev_colors)
          } else {
            plot_list[[m_var]] <- plot_list[[m_var]] +
          scale_fill_gradient(low = "#F0F0F0", high = "black",
                              breaks = c(1, ticks, max_val),
                              labels = c(1, ticks, max_val)) +
              guides(fill = guide_colourbar(reverse = TRUE))
            }
    } else if (tl %in% c("Core_or_not", "Host_group") & m_var == "Sample_ID") {
      country_plots <- list()
      for (country in unique(tible_list[[m_var]]$Country)) {
        country_plots[[country]] <- tible_list[[m_var]] %>%
          filter(Country == country) %>%
          ggplot(aes(x=.data[[m_var]], y=mean_tpm, fill = as.factor(group))) +
          geom_col() +
          ggtitle(paste0(country, " - ", title_prefix, hg_or_core, ": \"", hg,"\"")) +
          labs(fill=tl) +
          theme(axis.text.x = element_text(angle = 45, hjust=1)) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()) +
          scale_x_discrete(expand = c(0.025, 0)) +
          geom_text(aes(label=in_group, y = 1.01, vjust = 0, angle = 45)) +
          facet_wrap(Season~Hive_ID, scales = "free_x", nrow = 1)
        
        if (hg_or_core == "Host_group") {
          country_plots[[country]] <- country_plots[[country]] +
            scale_fill_manual(values = host_colors)
        }
      }
      plot_list[[m_var]] <- wrap_plots(country_plots, ncol = 1)
    } else {
        plot_list[[m_var]] <- tible_list[[m_var]] %>%
          ggplot(aes(x=.data[[m_var]], y=mean_tpm, fill = as.factor(group))) +
          geom_col() +
          ggtitle(paste0(title_prefix, hg_or_core, ": \"", hg,"\"")) +
          labs(fill=tl) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())
      }
    if (hg_or_core == "Host_group") {
      plot_list[[m_var]] <- plot_list[[m_var]] +
        scale_fill_manual(values = host_colors)
    }
    # if (m_var=="Sample_ID") {
    #   plot_list[[m_var]] <- plot_list[[m_var]] +
    #     theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    #     scale_x_discrete(expand = c(0.025, 0)) +
    #     geom_text(aes(label=in_group, y = 1.01, vjust = 0, angle = 45))
    # }
  }
  return(list(plots = plot_list, tibbles = tible_list))
}

prevalence_histogram <- function(abtable, plot_title) {
  plot_tbl <- abtable %>%
    pivot_longer(-contig) %>%
    mutate(value = ifelse(value>0, 1, 0)) %>%
    group_by(contig) %>%
    summarise(prevalence_abs = sum(value),
              prevalence_prop = sum(value) / (ncol(abtable)-1))
  hist_plot <- plot_tbl %>%
    ggplot(aes(x=prevalence_abs)) +
    geom_histogram(binwidth = 1, color="white") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    scale_x_continuous(breaks = 1:(ncol(abtable)-1)) +
    labs(title=plot_title)
  return(list(table = plot_tbl, plot = hist_plot))
}

prevalence_bar_plot <- function(abtable, tl, hg, meta_vars, title_prefix="", threshold_for_other=1) {
  label_for_other <- paste0("other (<", threshold_for_other,")")
  tax <- colnames(abtable)[1]
  title_suffix=paste0(" (prevalence >= ",threshold_for_other,")")
  plot_list <- list()
  for (m_var in meta_vars) {
    metadata_filt <- metadata %>%
      select(Sample_ID, all_of(m_var))
    
    plot_tbl <- abtable %>%
    rename(group = all_of(tax)) %>%
    pivot_longer(-group, names_to = "Sample_ID") %>%
    mutate(value = ifelse(value>0, 1, 0))  %>% 
    left_join(., metadata_filt, by="Sample_ID") %>%
    group_by(.data[[m_var]], group) %>%
    summarise(prevalence = sum(value), .groups="drop") %>%
    group_by(group) %>%
    mutate(group =
             ifelse(sum(prevalence) < threshold_for_other,
                    label_for_other,
                    group)) %>% 
    group_by(.data[[m_var]], group) %>%
    mutate(prevalence = sum(prevalence)) %>%  # This part is only to
    ungroup() %>%                         # merge identically-
    distinct()                       # colored shapes.
    
    # To plot in order of prevalence
    prev_oder <- plot_tbl %>%
      group_by(group) %>%
      summarise(total=sum(prevalence)) %>%
      arrange(desc(total)) %>%
      select(group) %>%
      unlist(use.names = FALSE)
    
    plot_list[[m_var]] <- plot_tbl %>%
      mutate(group = factor(group, levels=prev_oder)) %>%
      ggplot(aes(x=group, fill=.data[[m_var]], y=prevalence)) +
      geom_col() +
      ggtitle(paste0(title_prefix, tl, " - Host group: \"", hg,"\"", title_suffix)) +
      labs(fill=m_var) +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            plot.margin = margin(10, 10, 10, 100)) +
      scale_fill_manual(values = color_vector) #+
      # scale_x_discrete(expand = c(0.5, 0))  # Adjust the first value as needed
  }
  wrap <- wrap_plots(plot_list, ncol=1)
  return(wrap)
}

upset_country <- function(abtable) {
  tax <- colnames(abtable)[1]
  sets <- abtable %>%
    rename(group = all_of(tax)) %>%
    pivot_longer(-group) %>%
    filter(value > 0) %>%
    select(-value) %>%
    group_by(name) %>%
    summarise(ids = list(group)) %>%
    deframe()
  
  # m <- make_comb_mat(sets)
  # degrees <- comb_degree(m) %>%
  #   as.vector() %>%
  #   table()
  
  colors <- brewer.pal(8, "Dark2")
  rep_colors <- unlist(mapply(rep, colors, c(1,7,20,31,42,38,27,8)), use.names = FALSE) 
  
  p <- upset(fromList(sets), order.by = "degree", nsets= 8, nintersects = NA, main.bar.color = rep_colors,
             mainbar.y.label = "Number of phages")
  return(p)
}
