library(tidyverse)

completeness_cutoff <- 40

classification <- readRDS("output/R/R_variables/classification.RDS")


prediction_result <- read.delim("~/Library/CloudStorage/OneDrive-KULeuven/PhD/Virome/BPhage/output/lifestyle/replidec/BC_predict.summary") %>%
  tibble() %>%
  rename(contig = sample_name) %>%
  select(contig, final_label) %>%
  filter(!str_detect(contig, "Blank")) %>%
  mutate(original_contig = str_extract(contig, "^([^_]+_){10}[^_]+")) %>%
  left_join(., classification, by = join_by(original_contig == contig)) %>%
  select(contig, original_contig, final_label, Class, Core)

# prediction_result <- read.delim("output/lifestyle/bacphlip/bphage.fasta.bacphlip") %>%
# prediction_result <- read.delim("output/lifestyle/bacphlip/bphage_and_extended.fasta.bacphlip") %>%
#   tibble() %>%
#   rename("contig" = "X") %>%
#   filter(!str_detect(contig, "Blank")) %>%
#   mutate(original_contig = str_extract(contig, "^([^_]+_){10}[^_]+")) %>% #,
#          # extension = str_remove(contig, "^([^_]+_){10}[^_]+"))
#   left_join(., classification, by = join_by(original_contig == contig)) %>%
#   select(contig, original_contig, Virulent, Temperate, Class, Core)

extended <- read.delim("output/core_contig_refinement/extended_contigs_checkv/quality_summary.tsv") %>%
  tibble() %>%
  rename(contig = contig_id) %>%
  left_join(., prediction_result, by = "contig") %>%
  # select(contig, original_contig, completeness, Virulent, Temperate, Class, Core)
  select(contig, original_contig, completeness, final_label, Class, Core)

unextended <- classification %>%
  tibble() %>%
  filter(!contig %in% extended$original_contig) %>%
  select(contig, completeness) %>%
  left_join(., prediction_result, by = "contig") %>%
  # select(contig, original_contig, completeness, Virulent, Temperate, Class, Core)
  select(contig, original_contig, completeness, final_label, Class, Core)

joined_tibble_all <- rbind(unextended, extended) %>%
  filter(!is.na(Core)) %>%
  # mutate(Lifestyle = ifelse(Virulent > Temperate, "Lytic", "Temperate")) %>%
  filter(completeness >= completeness_cutoff)

lifestyle_counts_caudos <- joined_tibble_all %>% 
  filter(Class == "Caudoviricetes") %>%
  # group_by(Lifestyle, Core) %>%
  group_by(final_label, Core) %>%
  summarise(count = n(), .groups = "drop") %>%
  # mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
  mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
  arrange(Core) %>%
  mutate(subset = "Caudoviricetes")
  
lifestyle_tibble <- joined_tibble_all %>%
  # group_by(Lifestyle, Core) %>%
  group_by(final_label, Core) %>%
  summarise(count = n(), .groups = "drop") %>%
  # mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
  mutate(final_label = factor(final_label, levels = c("Virulent", "Chronic", "Temperate"))) %>%
  arrange(Core) %>%
  mutate(subset = "all_phages") %>%
  rbind(., lifestyle_counts_caudos)

lifestyle_bars <- list()
# lifestyle_colors <- c("#1C3A3A", "#FFC300")
lifestyle_colors <- c("#1C3A3A", "#8B4513", "#FFC300")
for (set in c("all_phages", "Caudoviricetes")) {
  lifestyle_bars[[set]] <- lifestyle_tibble %>%
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

system("mkdir -p output/R/lifestyle/")
for (set in names(lifestyle_bars)) {
  # ggsave(paste0("output/R/lifestyle/", set, ".completeness_", completeness_cutoff,".pdf"),
  ggsave(paste0("output/R/lifestyle/replidec_", set, ".completeness_", completeness_cutoff,".pdf"),
         lifestyle_bars[[set]], width = 6, height = 5)
}
write_csv(lifestyle_tibble,
          # paste0("output/R/lifestyle/lifestyle_stats.completeness_", completeness_cutoff, ".csv"))
          paste0("output/R/lifestyle/replidec_lifestyle_stats.completeness_", completeness_cutoff, ".csv"))
