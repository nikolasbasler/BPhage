library(vegan)
library(tidyverse)

metadata <- readRDS("data/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()
present_in_all_countries <- read_lines("data/core_contigs.txt")

phages <- list()
phages$core <- present_in_all_countries
phages$all <- phage_tpm$contig


## Individual samples and different gut parts 
phage_presence_mat <- list()
gp_list <- list(mid = "mid", ile = "ile", rec = "rec", all = c("mid", "ile", "rec"))
for (set in names(phages)) {
  for (gut_set in names(gp_list)) {
    
    phage_presence_mat[[set]][[gut_set]] <- phage_tpm %>%
      pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
      left_join(., metadata[c("Sample_ID", "Gut_part")], by = "Sample_ID") %>%
      filter(
        Gut_part %in% gp_list[[gut_set]],
        contig %in% phages[[set]],
      ) %>%
      group_by(contig, Sample_ID) %>%
      summarise(phage_presence = ifelse(sum(tpm) > 0, TRUE, FALSE), .groups = "drop") %>%
      pivot_wider(names_from = contig, values_from = phage_presence, values_fill = FALSE) %>%
      column_to_rownames("Sample_ID") %>%
      as.matrix()
  }
}

## Bee pools
for (set in names(phages)) {
  phage_presence_mat[[set]]$Bee_pool <- phage_tpm %>%
    pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
    left_join(., metadata[c("Sample_ID", "Bee_pool")], by = "Sample_ID") %>%
    filter(
      contig %in% phages[[set]],
    ) %>%
    group_by(contig, Bee_pool) %>%
    summarise(phage_presence = ifelse(sum(tpm) > 0, TRUE, FALSE), .groups = "drop") %>%
    pivot_wider(names_from = contig, values_from = phage_presence, values_fill = FALSE) %>%
    column_to_rownames("Bee_pool") %>%
    as.matrix()
  
}

## Accumulation curves
accumulation_plot <- list()
total_observed <- tibble()
accum_tibble <- list()
for (set in names(phage_presence_mat)) {
  for (met_v in names(phage_presence_mat[[set]])) {
  
    accum <- specaccum(phage_presence_mat[[set]][[met_v]],
                       method       = "random",
                       permutations = 100)
    
    accum_tibble[[set]][[met_v]] <- tibble(sites = accum$sites,
                           richness = accum$richness,
                           sd = accum$sd
    ) %>%
      mutate(lower = richness - sd,
             upper = richness + sd)

    total_observed <- tibble(slice = paste0(set, "_", met_v),
           total_observed = max(accum_tibble[[set]][[met_v]]$richness)) %>%
      rbind(total_observed, .)
    
    accumulation_plot[[set]][[met_v]] <- ggplot(accum_tibble[[set]][[met_v]], aes(x = sites, y = richness)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_point(size = 0.5) +
      labs(
        x     = "Number of samples",
        y     = "Accumulated number of unique phages" ) +
      theme_minimal(base_size = 14) +
      ggtitle(paste0(set, " phages - gut part: ", met_v))
    
  }
}

#### Save files
system("mkdir -p output/R/accumulation_curves/")
for (set in names(accumulation_plot)) {
  for (met_v in names(accumulation_plot[[set]])) {
    ggsave(paste0("output/R/accumulation_curves/accumulation_curve.", set, ".", met_v,".pdf"),
           accumulation_plot[[set]][[met_v]],
           width = 6, height = 6)
    
    write_tsv(
      accum_tibble[[set]][[met_v]],
      paste0("output/R/accumulation_curves/accumulation_curve.", set, ".", met_v,".tsv")
      )
  }
}
write_csv(total_observed, "output/R/accumulation_curves/total_observed.csv")
