library(ggVennDiagram)
library(patchwork)
library(tidyverse)

present_in_all_countries <- read_lines("data/core_contigs.txt")
# classification_core <- readRDS("output/R/R_variables/classification.RDS") %>%
#   filter(Core == "yes")


# LOAD THE OTHER DATASET
bphage_and_others_clusters <- read.delim("output/bphage_and_others_clusters.tsv", header=FALSE) %>%
  rename(representative = V1,
         member = V2)
  # filter(str_detect(member, "NODE") & 
  #          (str_detect(member, "Deboutte") | str_detect(member, "Busby") | str_detect(member, "Bonilla")
  #           )) %>%

clusters <- list()

clusters$all_BPhage <- bphage_and_others_clusters %>%
  separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")

# This is still not correct. It excludes all genomes that are exclusive to one of the
# other datasets!!
clusters$core_BPhage <- bphage_and_others_clusters %>%
  filter(str_detect(member, str_c(present_in_all_countries, collapse = "|"))) %>%
  separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")

# presence_in_datasets <- list()
venn <- list()
for (core_or_not in names(clusters)) {
  presence_in_datasets <- clusters[[core_or_not]] %>%
    mutate(across(
      starts_with("member_"),
      ~ case_when(
        str_detect(., "NODE") ~ "BPhage",
        str_detect(., "Deboutte") ~ "Deboutte",
        str_detect(., "Busby") ~ "Busby",
        str_detect(., "Bonilla") ~ "Bonilla",
        TRUE ~ .
      ))) %>%
    pivot_longer(-representative, values_drop_na = TRUE, names_to = "member", values_to = "dataset") %>%
    select(-member) %>%
    distinct()
  
  genomes_in_dataset <- list()
  for (set in c("BPhage", "Deboutte", "Busby", "Bonilla")) {
    genomes_in_dataset[[set]] <- presence_in_datasets %>%
      filter(dataset == set) %>%
      select(representative) %>%
      unlist(use.names = FALSE)
  }
  
 venn[[core_or_not]] <- ggVennDiagram(genomes_in_dataset) +
    theme(legend.position = "none")
}

