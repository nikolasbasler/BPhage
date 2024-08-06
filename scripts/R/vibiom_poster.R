library(tidyverse)
library(ggpubr)
library(patchwork)
library(stringr)

bee_brown <- "#5E280C"

metadata <- read.csv("data/metadata.csv") %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  mutate(Season = ifelse(Season == "spr", "spring", Season),
         Season = ifelse(Season == "sum", "summer", Season),
         Season = ifelse(Season == "aut", "autumn", Season)) %>%
  mutate(Season = factor(Season, levels = c("spring", "summer", "autumn"))) %>%
  mutate(Gut_part = ifelse(Gut_part == "mid", "midgut", Gut_part),
         Gut_part = ifelse(Gut_part == "ile", "ileum", Gut_part),
         Gut_part = ifelse(Gut_part == "rec", "rectum", Gut_part)) %>%
  mutate(Gut_part = factor(Gut_part, levels = c("midgut", "ileum", "rectum"))) %>%
  mutate(Health = factor(Health, levels = c("very_bad", "bad", "good", "very_good")))
sample_order <- metadata %>% 
  arrange(Country, Season, Hive_ID, Gut_part) %>%
  select(Sample_ID) %>%
  unlist(use.names = FALSE)
metadata <- metadata %>%
  mutate(Sample_ID = factor(Sample_ID, levels=sample_order))
row.names(metadata) <- metadata$Sample_ID


present_in_all_countries <- read_lines("data/core_contigs.txt")


####################  ALPHA ##################################
alpha_stats <- read.csv("output/R/alpha/alpha.contig.csv")
poster_alpha <- list()
for (metavar in c("Country", "Season", "Gut_part")) {
  ylabel = NULL
  
  if (metavar == "Country") {
    p_val = "7.9e-17"
    ylabel = "Hill-Shannon"
  }
  if (metavar == "Season") {
    p_val <- "0.6"
  }
  if (metavar == "Gut_part") {
    p_val <- "3.4e-05"
  }

  test_text <- paste0("p = ", p_val)
  
  poster_alpha[[metavar]] <- left_join(alpha_stats, metadata, by = "Sample_ID") %>%
    select(Sample_ID, Hill_Shannon, all_of(metavar)) %>%
    ggplot(aes(x = .data[[metavar]], y = Hill_Shannon)) +
    geom_boxplot(fill = "#ef8f01") +
    ggtitle(test_text) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = str_replace(metavar, "_", " ") , y = ylabel) +
    coord_cartesian(ylim = c(4.3, 91)) +
    theme(axis.ticks.y=element_blank(),
          panel.grid.minor = element_blank())
  
  if (metavar != "Country") {
    poster_alpha[[metavar]] <- poster_alpha[[metavar]] +
      theme(axis.text.y = element_blank())
  }
}

alpha_plot <- poster_alpha$Country + poster_alpha$Gut_part + poster_alpha$Season + 
  plot_layout(widths = c(1.3, 1, 1), ncol=3)
alpha_plot
# 
# ggsave(paste0("output/R/vibiom_poster/alpha.pdf"),
#        alpha_plot, width=6, height=3)

####################  PREVALENCE ##################################

prevalence_stats <- read.csv("output/R/prevalence/prevalence.Countries.csv")

color_vector <- "#ef8f01"

prev_plot <- prevalence_stats %>%
  group_by(prevalence_abs) %>%
  summarise(count = n()) %>%
  mutate(prevalence_abs = factor(prevalence_abs, levels = unique(prevalence_abs))) %>%
  ggplot(aes(x = prevalence_abs, y = count)) +
  geom_col(color = "black", fill = color_vector) +
  geom_text(aes(label=count), vjust=-0.5) +
  labs(x = "Number of countries", y = "Phage count") +
  coord_cartesian(ylim = c(60, 1525)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor = element_blank()) 
prev_plot
# ggsave(paste0("output/R/vibiom_poster/prevalence.pdf"),
#        prev_plot, width = 3, height = 2.4)

 ####################  MEAN REL ABUND ##################################

average_TPM_stats <- read.csv("output/R/relative_abundance_overall/average.TPM.Host_groups.Country.csv")
color_vector <- c("black","#ef8f01")

overall_mean <- 0.196

averate_TPM_plot <- average_TPM_stats %>%
  mutate(Country = factor(Country, levels = levels(metadata$Country))) %>%
  mutate(group = ifelse(group == "all","non-core", group),
         group = ifelse(group == "core", "core   ", group)) %>%
  ggplot(aes(x = Country, y = mean_tpm, fill = factor(group, levels = c("non-core", "core   ")))) +
  geom_col(color = "black") +
  scale_fill_manual(values = color_vector) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Mean relative abundance", fill = NULL) +
  coord_cartesian(ylim = c(0.043, 0.955)) +
  theme(panel.background = element_blank(),
        legend.text = element_text(angle = 90, vjust=0.5),
        legend.margin=margin(0,-4,0,-12),
        axis.text.y = element_text(margin = margin(r = -4)),
        axis.title.y = element_text(vjust=0),
        axis.ticks.y = element_blank()
        ) +
  guides(fill = guide_legend(label.position = "top", keywidth = 1, keyheight = 1))

averate_TPM_plot

# ggsave(paste0("output/R/vibiom_poster/TPM.pdf"),
#        averate_TPM_plot, width = 3, height = 2.25)

####################  HOSTS ##################################

iphop_prediction_genus <- read.csv("output/host_prediction/iphop_output_core/Host_prediction_to_genus_m75.csv")
confident_host <- iphop_prediction_genus %>%
  filter(Confidence.score >= 90) %>% 
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score)) %>%
  select(Virus, Host.genus) %>%
  mutate(only_genus = str_split_i(Host.genus, ';', 6)) %>%
  mutate(only_genus = str_replace_all(only_genus, "g__", "")) %>%
  mutate(only_genus = ifelse (only_genus == "Apilactobacillus", "other", only_genus),
         only_genus = ifelse (only_genus == "Acinetobacter", "other", only_genus),
         only_genus = ifelse (only_genus == "Frischella", "other", only_genus)) %>%
  select(Virus, only_genus)

unknown_host <- confident_host %>%
  select(Virus) %>%
  setdiff(tibble(Virus=present_in_all_countries), .)

host_tibble <- bind_rows(confident_host, tibble(unknown_host, only_genus = "unkown"))

host_pie_colors <- rev(c("#FFDAB9", "#FFA07A", "#DAA520", "#D2691E", "#8B4513", "#555555", "lightgrey"))

host_pie <- host_tibble %>%
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

# ggsave(paste0("output/R/vibiom_poster/host.pdf"),
#        host_pie, width = 4.5, height = 4.5)

####################  LIFESTYLE ##################################

lifestyle_tibble <- read.delim("output/lifestyle/bacphlip/bphage.fasta.bacphlip") %>%
  rename("contig" = "X") %>%
  filter(contig %in% present_in_all_countries)

lifestyle_pie_colors <- c("black", "#ef8f01")

lifestyle_bar <- lifestyle_tibble %>%
  mutate(Lifestyle = ifelse(Virulent > Temperate, "Lytic", "Temperate")) %>%
  group_by(Lifestyle) %>%
  summarise(count = n()) %>%
  mutate(Lifestyle = factor(Lifestyle, levels = c("Lytic", "Temperate"))) %>%
  ggplot(aes(x = "", y = count, fill = Lifestyle)) +
  geom_bar(stat = "identity", color= "black") +
  theme_void() +
  theme(legend.text.align = 1,
        legend.margin=margin(0,20,0,-20)) +
  scale_fill_manual(values = lifestyle_pie_colors) +
  geom_text(aes(label = c("46", "51")), color = c("#ef8f01", "black"), size = 16,
            position = position_stack(vjust = 0.5))

# ggsave(paste0("output/R/vibiom_poster/lifestyle.pdf"),
#        lifestyle_bar, width = 3, height = 10)

####################  AMGs ##################################

pharokka_final_merged <- read.delim("output/annotation/old_pharokka/core_contigs/core_contigs.fasta_cds_final_merged_output.tsv")

moron_bar_colors <- rev(c("lightgrey", "#ef8f01", "#555555"))

phrog_bar <- pharokka_final_merged %>%
  group_by(category) %>%
  summarise(count = n()) %>%
  mutate(collapsed_cat = "Other PHROG\ncategories") %>%
  mutate(collapsed_cat = ifelse(category == "unknown function", "unknown function", collapsed_cat),
         collapsed_cat = ifelse(category == "moron, auxiliary metabolic gene and host takeover", '"moron", AMG\nand host takeover', collapsed_cat)
         ) %>%
  group_by(collapsed_cat) %>%
  summarise(coll_count = sum(count)) %>%
  mutate(collapsed_cat = factor(collapsed_cat, levels = c("Other PHROG\ncategories", '"moron", AMG\nand host takeover', "unknown function"))) %>%
  ggplot(aes(x = "", y = coll_count, fill = collapsed_cat)) +
  geom_bar(stat = "identity") +
  theme_void() +
  theme(legend.text.align = 1,
        legend.margin=margin(0,20,10,0)) +
  scale_fill_manual(values = moron_bar_colors) +
  labs(fill = "PHROG category") +
  theme(legend.spacing.x = unit(0.25, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))

phrog_bar 

# ggsave(paste0("output/R/vibiom_poster/phrog.pdf"),
#        phrog_bar, width = 2.5, height = 6)


moron_tibble <- pharokka_final_merged %>%
  filter(category == "moron, auxiliary metabolic gene and host takeover") %>% 
  group_by(annot) %>%
  summarise(count = n()) %>%
  mutate(group = NA) %>%
  mutate(group = ifelse(annot == "membrane protein" |
                          annot == "membrane associated protein",
                        "Membrane protein", group),
         group = ifelse(annot == "phosphoadenosine phosphosulfate reductase", "Sulfur metabolism", group),
         group = ifelse(annot == "gam-like host nuclease inhibitor" |
                          annot == "anti-restriction protein" |
                          annot == "anti-sigma factor",
                        "Host takeover", group),
         group = ifelse(annot == "Doc-like toxin"|
                          annot == "MazF-like growth inhibitor" |
                          annot == "toxin" |
                          annot == "toxin-antitoxin system HicB-like" |
                          annot == "Lar-like restriction alleviation protein" |
                          annot == "RelE-like toxin" |
                          annot == "plasmid antitoxin with HTH domain" |
                          annot == "ribonuclease toxin of AT system",
                        "Toxin-antitoxin system", group),
         group = ifelse(annot == "abortive infection resistance protein" |
                          annot == "superinfection exclusion" |
                          annot == "SieB superinfection exclusion",
                        "Superinfection exclusion", group),
         group = ifelse(annot == "beta-lactamase-inhibitor protein BLIP", 
                        "Antibiotic resistance", group),
         group = ifelse(annot == "queuine tRNA-ribosyltransferase" |
                          annot == "PnuC-like nicotinamide mononucleotide transport",
                        "other", group)
         ) %>%
  group_by(group) %>%
  summarise(group_count = sum(count)) %>%
  mutate(group = factor(group, levels = rev(c("Membrane protein", "Toxin-antitoxin system",
                                                      "Sulfur metabolism", "Host takeover",
                                                      "Superinfection exclusion", "Antibiotic resistance",
                                                      "other"))))

sum(moron_tibble$group_count)

pharokka_final_merged %>%
  filter(category == "moron, auxiliary metabolic gene and host takeover") %>% 
  select(contig) %>%
  distinct() %>%
  nrow()

moron_pie_colors <- c("#555555", "#FFDAB9", "#FFA07A", "#DAA520", "#D2691E", "#8B4513", "#5E280C")
moron_pie <- 
  moron_tibble %>%
  ggplot(aes(x = "", y = group_count, fill = group)) +
  geom_bar(stat = "identity", color= "black") +
  coord_polar("y") +
  theme_void() +
  guides(fill = guide_legend(reverse=TRUE)) +
  labs(fill = 'Protein group') +
  theme(legend.margin=margin(0,2,0,-20)) +
  scale_fill_manual(values = moron_pie_colors)

moron_pie 

ggsave(paste0("output/R/vibiom_poster/moron.pdf"),
       moron_pie, width = 4.5, height = 4.5)
  

#### membrane protein:
#   phrog_1265: ???
#   phrog_1486: ???
#   phrog_3996: ???
#   phrog_502: ???
#   phrog_7180: tellurite resistance protein TerC (KO)

#### phosphoadenosine phosphosulfate reductase: Sulfur metabolism (KO)
#### gam-like host nuclease inhibitor: Host take-over ? Nuclease inhibition (Pfam)
#### Doc-like toxin: Toxin-antitoxin system (KO)
#### MazF-like growth inhibitor: Toxin-antitoxin system (KO)
#### abortive infection resistance protein: Other (Superinfection exclusion? Pfam)
#### beta-lactamase-inhibitor protein BLIP: Antibiotic resistence (Pfam)
#### membrane associated protein: ???
#### queuine tRNA-ribosyltransferase: ???
#### superinfection exclusion: Superinfection exclusion (no KO)
#### toxin: Toxin-antitoxin system (KO)
#### toxin-antitoxin system HicB-like: Toxin-antitoxin system (KO)
#### Lar-like restriction alleviation protein: Toxin-antitoxin system (KO)
#### PnuC-like nicotinamide mononucleotide transport: Other
#### RelE-like toxin: Toxin-antitoxin system (KO)
#### SieB superinfection exclusion: Superinfection exclusion
#### anti-restriction protein: Host takeover?
#### anti-sigma factor: Host takeover?
#### plasmid antitoxin with HTH domain: Toxin-antitoxin system (PHROG)
#### ribonuclease toxin of AT system: Toxin-antitoxin system (PHROG)


  
  
  
  
  
  
  
  
  
  
  
  
#










