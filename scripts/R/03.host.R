library(tidyverse)

# color_vector <- c("#FFDAB9", "#FFA07A", "#FFC300","#ef8f01", "#D2691E", "#8B4513", "#1C3A3A", "black", "#555555", "lightgrey")

host_pie_colors <- c("Bifidobacterium" = "#FFDAB9",
                     "Lactobacillus" = "#FFA07A",
                     "Snodgrassella" = "#FFC300",
                     "Bombilactobacillus" = "#ef8f01",
                     "Gilliamella" = "#8B4513",
                     "Frischella" = "#336B6B",
                     "Commensalibacter" = "#285353",
                     "Bartonella" = "#1C3A3A",
                     "Bombella" = "#112222",
                     "other" = "#555555",
                     "unknown" = "lightgrey")

classification <- readRDS("output/R/R_variables/classification.RDS")
all_contigs <- classification$contig %>% as.character()
present_in_all_countries <- read_lines("data/core_contigs.txt")

iphop_prediction <- list()
# iphop_prediction$all <- read.csv("output/host_prediction/iphop_bphage_and_others_custom_database/iphop_v1.3.3_db_Aug_2023_Host_prediction_to_genome_m90.csv") %>%
#   filter(str_starts(Virus, "NODE")) %>%
#   filter(!str_detect(Virus, "Blank_pool")) # This one should have been filtered out earlier
# iphop_prediction$all <- read.csv("output/host_prediction/iphop_bphage_and_others_custom_database/iphop_v1.3.3_db_Aug_2023_Host_prediction_to_genus_m90.csv") %>%
#   filter(str_starts(Virus, "NODE")) %>%
#   filter(!str_detect(Virus, "Blank_pool"))
# iphop_prediction$all <- read.csv("output/host_prediction/iphop_output_bphage/Host_prediction_to_genome_m75.csv") %>%
#   filter(!str_detect(Virus, "Blank_pool"))
iphop_prediction$all <- read.csv("output/host_prediction/iphop_output_bphage/Host_prediction_to_genus_m75.csv") %>% 
  filter(!str_detect(Virus, "Blank_pool"))

iphop_prediction$core <- iphop_prediction$all %>%
  filter(Virus %in% present_in_all_countries)
iphop_prediction$noncore <- iphop_prediction$all %>%
  filter(!Virus %in% present_in_all_countries)

confident_host <- list()
confident_host_genus <- list()
confident_host_genome <- list()
for (set in names(iphop_prediction)) {
  confident_host[[set]] <- iphop_prediction[[set]] %>%
    filter(Confidence.score >= 90) %>% 
    group_by(Virus) %>%
    filter(Confidence.score == max(Confidence.score)) %>%
    ungroup()
  
  confident_host_genus[[set]] <- confident_host[[set]] %>%
    group_by(Virus) %>%
    slice(1) %>% # In case several genomes have the same confidence score, only keep the first one.
    ungroup() %>%
    
    select(Virus, Host.genus) %>%
    mutate(Genus = sapply(str_split(Host.genus, ";"), function(x) {
      genus_field <- x[grepl("^g__", x)]
      if(length(genus_field) > 0) return(genus_field) else return("g__")
    })) %>%
    # select(Virus, Host.taxonomy) %>%
    # mutate(Genus = sapply(str_split(Host.taxonomy, ";"), function(x) {
    #   genus_field <- x[grepl("^g__", x)]
    #   if(length(genus_field) > 0) return(genus_field) else return("g__")
    # })) %>% 
    filter(Genus != "g__") %>%
    mutate(Genus = str_replace_all(Genus, "g__", "")) %>%
    mutate(Genus = str_replace_all(Genus, "_.$", "")) %>% 
    select(Virus, Genus)
  
  confident_host_genus[[set]]
  
  # confident_host_genome[[set]] <- confident_host[[set]] %>%
  #   group_by(Virus) %>%
  #   filter(if(any(grepl("_bin.", Host.genome))) grepl("_bin.", Host.genome) else TRUE) %>% # If one of bacterial BPhage MAGs is among the hosts, keep that instead of the standard database hit
  #   slice(1) %>%
  #   ungroup()
}

unknown_host <- list()
unknown_host$all <- setdiff(all_contigs, confident_host_genus$all$Virus) %>%
  tibble(Virus= .)
unknown_host$core <- confident_host_genus$core %>%
  select(Virus) %>%
  setdiff(tibble(Virus=present_in_all_countries), .)
unknown_host$noncore <- classification %>%
  filter(Core =="no",
         !contig %in% confident_host_genus$noncore$Virus) %>%
  select(contig) %>%
  rename(Virus = contig) %>%
  tibble()

host_pie <- list()
host_tibble <- list()
host_group <- list()
all_hosts <- list()
for (set in c("all", "core", "noncore")) {
  all_hosts[[set]] <- bind_rows(confident_host_genus[[set]], tibble(unknown_host[[set]], Genus = "unknown"))
  host_group[[set]] <- all_hosts[[set]] %>%
    mutate(Genus = ifelse(Genus == "", "unknown", Genus)) %>%
    mutate(Genus = ifelse(!Genus %in% c("Gilliamella", "Lactobacillus",
                                        "Bifidobacterium", "Snodgrassella",
                                        "Bombilactobacillus", "Bartonella", 
                                        "Frischella", "Bombella" , "Commensalibacter", "unknown"), 
                          "other", Genus))
  
  host_tibble[[set]] <- host_group[[set]] %>%
    group_by(Genus) %>%
    summarise(count = n()) %>% 
    arrange(count) %>%
    ungroup()
  
  genus_order <- host_tibble[[set]] %>%
    filter(!Genus %in% c("other", "unknown")) %>% 
    
    mutate(core_bacterium = ifelse(Genus %in% c("Gilliamella", "Lactobacillus", "Bifidobacterium", "Bombilactobacillus", "Snodgrassella"), TRUE, FALSE)) %>%
    arrange(core_bacterium) %>%
    
    
    select(Genus) %>%
    unlist(use.names = FALSE)

  host_pie[[set]] <- host_tibble[[set]] %>%
    mutate(Genus = factor(Genus, levels = c("unknown", "other", genus_order))) %>%
    ggplot(aes(x = "", y = count, fill = Genus)) +
    geom_bar(stat = "identity", color= "black") +
    coord_polar("y") +
    theme_void() +
    guides(fill = guide_legend(reverse=TRUE)) +
    labs(fill = "Host genus") +
    theme(legend.margin=margin(0,2,0,-20),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = host_pie_colors) +
    ggtitle(set)
}

system("mkdir -p output/R/host_pies")
for (pie in names(host_pie)) {
  ggsave(paste0("output/R/host_pies/hosts.", pie, ".pdf"),
         host_pie[[pie]], height = 6, width = 6)
  write_csv(host_tibble[[pie]],
            paste0("output/R/host_pies/hosts.", pie, ".csv"))
  write_csv(all_hosts[[pie]],
            paste0("output/R/host_pies/all_hosts.", pie, ".csv"))
}

# Written do data/ for convenience to avoid back tracking. So it can be used in the main analysis R script.
# host_group$all %>%
#   rename(contig = Virus) %>%
#   rename(Host_group = Genus) %>%
#   write_csv("data/host_groups.csv")
 