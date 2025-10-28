library(patchwork)
library(tidyverse)

lifestyle_colors <- c("Virulent" = "#1C3A3A", "Chronic" =  "#8B4513", "Temperate" ="#FFC300")

classification <- readRDS("data/classification.RDS")

replidec <- read.delim("output/lifestyle/replidec/BC_predict.summary") %>%
  tibble() 

extended_contigs <- read.delim("output/core_contig_refinement/extended_contigs_checkv/quality_summary.tsv") %>%
  tibble()

prediction_result <- replidec %>%
  rename(contig = sample_name) %>%
  select(contig, final_label) %>%
  filter(!str_detect(contig, "Blank")) %>%
  mutate(original_contig = str_extract(contig, "^([^_]+_){10}[^_]+")) %>%
  left_join(., classification, by = join_by(original_contig == contig)) %>%
  select(contig, original_contig, final_label, Class, Core)

extended <- extended_contigs %>%
  rename(contig = contig_id) %>%
  left_join(., prediction_result, by = "contig") %>% 
  select(contig, original_contig, completeness, final_label, Class, Core)

unextended <- classification %>%
  tibble() %>%
  filter(!contig %in% extended$original_contig) %>%
  select(contig, completeness) %>%
  left_join(., prediction_result, by = "contig") %>%
  select(contig, original_contig, completeness, final_label, Class, Core)

lifestyle_tibble <- list()
lifestyle_bars <- list()
lifestyle_bars_horizontal <- list()
cuts <- c("all" = 40, ">=75%" = 75, ">=90%" = 90)
for (cutoff in names(cuts)) {
  threshold <- cuts[[cutoff]]
  
  joined_tibble_all <- rbind(unextended, extended) %>%
    filter(!is.na(Core)) %>%
    filter(completeness >= threshold)
  
  lifestyle_counts_caudos <- joined_tibble_all %>% 
    filter(Class == "Caudoviricetes") %>%
    group_by(final_label, Core) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
    arrange(Core) %>%
    mutate(subset = "Caudoviricetes")
    
  lifestyle_tibble[[cutoff]] <- joined_tibble_all %>%
    group_by(final_label, Core) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
    arrange(Core) %>%
    mutate(subset = "all_phages") %>%
    rbind(., lifestyle_counts_caudos)
  
  for (set in c("all_phages", "Caudoviricetes")) {
    lifestyle_bars[[set]][[cutoff]] <- lifestyle_tibble[[cutoff]] %>%
      filter(subset == set) %>%
      group_by(Core) %>%
      mutate(proportion = count / sum(count)) %>%
      ggplot(aes(x = Core, y = proportion, fill = final_label)) +
      geom_bar(stat = "identity", color= "black") +
      scale_fill_manual(values = lifestyle_colors) +
      labs(title = paste0(set, " - completeness: ", cutoff),
           fill = "Lifestyle") +
      theme_minimal() +
      theme(legend.position = "bottom", legend.justification = "center") 
    lifestyle_bars_horizontal[[set]][[cutoff]] <- lifestyle_bars[[set]][[cutoff]] + 
      coord_flip() +
      scale_x_discrete(limits = rev)
  }
}

#####
# Save files

system("mkdir -p output/R/lifestyle")
for (cutoff in names(lifestyle_tibble)) {
  for_filename <- gsub("[>=%]", "", cutoff) # These indices were not my best idea...

  write_tsv(lifestyle_tibble[[cutoff]], paste0("output/R/lifestyle/replidec_completeness.", for_filename, ".tsv"))

  for (set in names(lifestyle_bars)) {
    ggsave(paste0("output/R/lifestyle/replidec.", set, ".", for_filename, ".vertical.pdf"),
           lifestyle_bars[[set]][[cutoff]], width = 3.75, height = 7)
    ggsave(paste0("output/R/lifestyle/replidec.", set, ".", for_filename, ".horizontal.pdf"),
           lifestyle_bars_horizontal[[set]][[cutoff]], width = 7, height = 3.75)
  }

}

# For convenience, to avoid backtracking
# all_predictions <- rbind(unextended, extended) %>%
#   select(original_contig, final_label) %>%
#   rename(contig = original_contig,
#          Lifestyle_replidec = final_label)
# 
# classification_strip <- classification %>%
#   select(-Lifestyle_replidec) # comment out if needed
# 
# new_classification_df <- left_join(classification_strip, all_predictions, by = "contig")
#
# write_csv(new_classification_df, "output/R/classification.csv")
