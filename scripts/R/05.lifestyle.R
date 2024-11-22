library(ggVennDiagram)
library(patchwork)
library(tidyverse)

# lifestyle_colors <- c("#1C3A3A", "#FFC300")
lifestyle_colors <- c("Virulent" = "#1C3A3A", "Chronic" =  "#8B4513", "Temperate" ="#FFC300")

classification <- readRDS("output/R/R_variables/classification.RDS")

cutoffs <- c(40, 75, 90)
for (completeness_cutoff in cutoffs) {
  
  prediction_result <- list()
  prediction_result$replidec <- read.delim("~/Library/CloudStorage/OneDrive-KULeuven/PhD/Virome/BPhage/output/lifestyle/replidec/BC_predict.summary") %>%
    tibble() %>%
    rename(contig = sample_name) %>%
    select(contig, final_label) %>%
    filter(!str_detect(contig, "Blank")) %>%
    mutate(original_contig = str_extract(contig, "^([^_]+_){10}[^_]+")) %>%
    left_join(., classification, by = join_by(original_contig == contig)) %>%
    select(contig, original_contig, final_label, Class, Core)
  
  # prediction_result$bacphlip <- read.delim("output/lifestyle/bacphlip/bphage.fasta.bacphlip") %>%
  prediction_result$bacphlip <- read.delim("output/lifestyle/bacphlip/bphage_and_extended.fasta.bacphlip") %>%
    tibble() %>%
    rename("contig" = "X") %>%
    filter(!str_detect(contig, "Blank")) %>%
    mutate(original_contig = str_extract(contig, "^([^_]+_){10}[^_]+")) %>% #,
           # extension = str_remove(contig, "^([^_]+_){10}[^_]+"))
    left_join(., classification, by = join_by(original_contig == contig)) %>%
    select(contig, original_contig, Virulent, Temperate, Class, Core) %>%
    mutate(final_label = ifelse(Virulent > Temperate, "Virulent", "Temperate"))
  
  lifestyle_tibble <- list()
  joined_tibble_all <- list()
  lifestyle_bars <- list()
  for (tool in names(prediction_result)) {
    extended <- read.delim("output/core_contig_refinement/extended_contigs_checkv/quality_summary.tsv") %>%
      tibble() %>%
      rename(contig = contig_id) %>%
      left_join(., prediction_result[[tool]], by = "contig") %>% 
      # select(contig, original_contig, completeness, Virulent, Temperate, Class, Core)
      select(contig, original_contig, completeness, final_label, Class, Core)
    
    unextended <- classification %>%
      tibble() %>%
      filter(!contig %in% extended$original_contig) %>%
      select(contig, completeness) %>%
      left_join(., prediction_result[[tool]], by = "contig") %>%
      # select(contig, original_contig, completeness, Virulent, Temperate, Class, Core)
      select(contig, original_contig, completeness, final_label, Class, Core)
    
    joined_tibble_all[[tool]] <- rbind(unextended, extended) %>%
      filter(!is.na(Core)) %>%
      # mutate(Lifestyle = ifelse(Virulent > Temperate, "Lytic", "Temperate")) %>%
      filter(completeness >= completeness_cutoff)
    
    lifestyle_counts_caudos <- joined_tibble_all[[tool]] %>% 
      filter(Class == "Caudoviricetes") %>%
      # group_by(Lifestyle, Core) %>%
      group_by(final_label, Core) %>%
      summarise(count = n(), .groups = "drop") %>%
      # mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
      mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
      arrange(Core) %>%
      mutate(subset = "Caudoviricetes")
      
    lifestyle_tibble[[tool]] <- joined_tibble_all[[tool]] %>%
      # group_by(Lifestyle, Core) %>%
      group_by(final_label, Core) %>%
      summarise(count = n(), .groups = "drop") %>%
      # mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
      mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
      arrange(Core) %>%
      mutate(subset = "all_phages") %>%
      rbind(., lifestyle_counts_caudos)
    
    for (set in c("all_phages", "Caudoviricetes")) {
      lifestyle_bars[[tool]][[set]] <- lifestyle_tibble[[tool]] %>%
        filter(subset == set) %>%
        group_by(Core) %>%
        mutate(proportion = count / sum(count)) %>%
        # ggplot(aes(x = Core, y = proportion, fill = Lifestyle)) +
        ggplot(aes(x = Core, y = proportion, fill = final_label)) +
        geom_bar(stat = "identity", color= "black") +
        scale_fill_manual(values = lifestyle_colors) +
        labs(title = paste0(set, " - completeness: >=", completeness_cutoff, "%"),
             fill = "Lifestyle")
    }
  }
  
  # Compare tools
  venn_by_style <- list()
  for (is_core in c("yes", "no")) {
    venn <- list()
    for (prediction in c("Virulent", "Temperate", "Chronic")) {
      contig_overlap <- list()
      for (tool in names(joined_tibble_all)) {
        contig_overlap[[is_core]][[tool]] <- joined_tibble_all[[tool]] %>%
          filter(final_label == prediction,
                 # Class == "Caudoviricetes",
                 Core == is_core) %>% 
          select(contig) %>%
          unlist(use.names = FALSE)
      }
      venn[[prediction]] <- ggVennDiagram(contig_overlap[[is_core]]) +
        theme(legend.position = "none") +
        coord_flip() +
        # ggtitle(paste0(completeness_cutoff, "% complete Caudos - core: ", is_core, " - ", prediction))
        ggtitle(paste0(prediction, " - ", completeness_cutoff, "% complete - core: ", is_core))
    }
    venn_by_style[[is_core]] <- wrap_plots(venn, ncol = 1)
  }
  
  system("mkdir -p output/R/lifestyle/")
  for (tool in names(lifestyle_bars)) {
    for (set in names(lifestyle_bars[[tool]])) {
      # ggsave(paste0("output/R/lifestyle/", set, ".completeness_", completeness_cutoff,".pdf"),
      ggsave(paste0("output/R/lifestyle/", tool, "_", set, ".completeness_", completeness_cutoff,".pdf"),
             lifestyle_bars[[tool]][[set]], width = 6, height = 5)
    }
    write_csv(lifestyle_tibble[[tool]],
              # paste0("output/R/lifestyle/lifestyle_stats.completeness_", completeness_cutoff, ".csv"))
              paste0("output/R/lifestyle/", tool, "_lifestyle_stats.completeness_", completeness_cutoff, ".csv"))
  }
  for (is_core in names(venn_by_style)) {
    ggsave(paste0("output/R/lifestyle/venn.core_", is_core, ".completeness_", completeness_cutoff, ".pdf"),
           venn_by_style[[is_core]], width = 6, height = 9)
  }
}
