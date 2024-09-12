library(tidyverse)

classification <- readRDS("output/R/R_variables/classification.RDS")

completeness_cutoff <- 75
joined_tibble_all <- read.delim("output/lifestyle/bacphlip/bphage.fasta.bacphlip") %>%
  rename("contig" = "X") %>%
  left_join(classification, ., by = "contig") %>%
  mutate(Lifestyle = ifelse(Virulent > Temperate, "Lytic", "Temperate")) %>%
  filter(completeness >= completeness_cutoff)

lifestyle_counts_caudos <- joined_tibble_all %>% 
  filter(Class == "Caudoviricetes") %>%
  group_by(Lifestyle, Core) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
  arrange(Core) %>%
  mutate(subset = "Caudoviricetes")
  
lifestyle_tibble <- joined_tibble_all %>%
  group_by(Lifestyle, Core) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
  arrange(Core) %>%
  mutate(subset = "all_phages") %>%
  rbind(., lifestyle_counts_caudos)

lifestyle_bars <- list()
lifestyle_colors <- c("black", "#ef8f01")
for (set in c("all_phages", "Caudoviricetes")) {
  lifestyle_bars[[set]] <- lifestyle_tibble %>%
    filter(subset == set) %>%
    group_by(Core) %>%
    mutate(proportion = count / sum(count)) %>%
    ggplot(aes(x = Core, y = proportion, fill = Lifestyle)) +
    geom_bar(stat = "identity", color= "black") +
    scale_fill_manual(values = lifestyle_colors) +
    labs(title = paste0(set, " - completeness: >", completeness_cutoff, "%"))
}

system("mkdir -p output/R/lifestyle/")
for (set in names(lifestyle_bars)) {
  ggsave(paste0("output/R/lifestyle/", set, ".completeness_", completeness_cutoff,".pdf"),
         lifestyle_bars[[set]], width = 6, height = 5)
}
write_csv(lifestyle_tibble, 
          paste0("output/R/lifestyle/lifestyle_stats.completeness_", completeness_cutoff, ".csv"))
