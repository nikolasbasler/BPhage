library(tidyverse)

classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")
  
iphop_prediction_genus <- list()
iphop_prediction_genus$core <- read.csv("output/host_prediction/iphop_output_core/Host_prediction_to_genus_m75.csv")
iphop_prediction_genus$noncore <- read.csv("output/host_prediction/iphop_output_bphage/Host_prediction_to_genus_m75.csv")

confident_host <- list()
for (set in names(iphop_prediction_genus)) {
  confident_host[[set]] <- iphop_prediction_genus[[set]] %>%
    filter(Confidence.score >= 90) %>% 
    group_by(Virus) %>%
    filter(Confidence.score == max(Confidence.score)) %>%
    select(Virus, Host.genus) %>%
    mutate(only_genus = str_split_i(Host.genus, ';', 6)) %>%
    mutate(only_genus = str_replace_all(only_genus, "g__", "")) %>%
    mutate(only_genus = ifelse (only_genus == "Apilactobacillus", "other", only_genus),
           only_genus = ifelse (only_genus == "Acinetobacter", "other", only_genus),
           only_genus = ifelse (only_genus == "Frischella", "other", only_genus)) %>%
    select(Virus, only_genus) %>%
    ungroup()
}
  
unknown_host <- list()
unknown_host$core <- confident_host$core %>%
  select(Virus) %>%
  setdiff(tibble(Virus=present_in_all_countries), .)
unknown_host$noncore <- classification %>%
  filter(Core =="no",
         !contig %in% confident_host$noncore$Virus) %>%
  select(contig) %>%
  rename(Virus = contig) %>%
  tibble()

host_pie <- list()
for (set in c("core", "noncore")) {
  host_tibble <- bind_rows(confident_host[[set]], tibble(unknown_host[[set]], only_genus = "unkown")) %>%
    mutate(only_genus = ifelse(only_genus == "", "unknown", only_genus)) %>%
    mutate(only_genus = ifelse(!only_genus %in% c("Gilliamella", "Lactobacillus",
                                                 "Bifidobacterium", "Snodgrassella",
                                                 "Bombilactobacillus", "unkown"), "other", only_genus))
  
  host_pie_colors <- rev(c("#FFDAB9", "#FFA07A", "#DAA520", "#D2691E", "#8B4513", "#555555", "lightgrey"))
  
  host_pie[[set]] <- host_tibble %>%
    group_by(only_genus) %>%
    summarise(count = n()) %>% 
    mutate(only_genus = factor(only_genus, levels = rev(c("Gilliamella", "Lactobacillus",
                                                          "Bifidobacterium", "Snodgrassella", 
                                                          "Bombilactobacillus", "other", "unkown")))) %>%
    ggplot(aes(x = "", y = count, fill = only_genus)) +
    geom_bar(stat = "identity", color= "black") +
    coord_polar("y") +
    theme_void() +
    guides(fill = guide_legend(reverse=TRUE)) +
    labs(fill = "Host genus") +
    theme(legend.margin=margin(0,2,0,-20)) +
    scale_fill_manual(values = host_pie_colors)
}

system("mkdir -p output/R/host_pies")
for (pie in names(host_pie)) {
  ggsave(paste0("output/R/host_pies/hosts.", pie, ".pdf"),
         host_pie[[pie]], height = 4, width = 4)
}
  