library(ggVennDiagram)
library(patchwork)
library(tidyverse)

present_in_all_countries <- read_lines("data/core_contigs.txt")

##### Contig overlap

clusters <- list()

clusters$all_BPhage <- read.delim("output/bphage_and_others_clusters.tsv", header=FALSE) %>%
  rename(representative = V1, member = V2) %>%
  separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")

clusters$core_BPhage <- read.delim("output/bphage_core_and_others_clusters.tsv", header=FALSE) %>%
  rename(representative = V1, member = V2) %>%
  separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")

conitg_overlap_venn <- list()
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
  
  conitg_overlap_venn[[core_or_not]] <- ggVennDiagram(genomes_in_dataset) +
    theme(legend.position = "none")
}

conitg_overlap_venn

##### Mapping

horizontal_coverage_threshold = 70
mean_depth_threshold = 1

SRA_to_study <- read.delim("data/other_datasets_SRA_accessions.tsv", header=FALSE) %>%
  rename(SRA = V1, study = V2)

stats.read.other_studies.reads <- read.delim("output/other_studies/stats.read.other_studies.reads.tsv") %>%
  select(SRA, Trimmed_pairs) %>%
  left_join(., SRA_to_study, by = "SRA") %>%
  group_by(study) %>%
  summarise(trimmed_read_pairs = sum(Trimmed_pairs))

read_counts <- read.delim("output/bphage_viper_output/read_stats.tsv") %>%
  mutate(Trimmed_pairs = Trimmed_R1_plus_R2 / 2) %>%
  mutate(study = "BPhage") %>%
  select(study, Trimmed_pairs) %>%
  group_by(study) %>%
  summarise(trimmed_read_pairs = sum(Trimmed_pairs)) %>%
  rbind(stats.read.other_studies.reads)

stats.other_studies.mapped_reads <- read.csv("output/other_studies/stats.other_studies.mapped_reads.csv") %>%
  pivot_longer(-contig, names_to = "SRA", values_to = "reads") 
stats.other_studies.horizontal_coverage <- read.csv("output/other_studies/stats.other_studies.horizontal_coverage.csv") %>%
  pivot_longer(-contig, names_to = "SRA" , values_to = "hzc")
stats.other_studies.mean_depth <- read.csv("output/other_studies/stats.other_studies.mean_depth.csv") %>%
  pivot_longer(-contig, names_to = "SRA", values_to = "depth")

filtered_ab_long <- inner_join(stats.other_studies.mapped_reads, stats.other_studies.horizontal_coverage, by = c("contig", "SRA")) %>%
  inner_join(., stats.other_studies.mean_depth, by = c("contig", "SRA")) %>%
  mutate(reads = ifelse(hzc < horizontal_coverage_threshold, 0, reads)) %>%
  mutate(reads = ifelse(depth < mean_depth_threshold, 0, reads)) %>%
  select(contig, SRA, reads)

presence_absence <- filtered_ab_long %>%
  filter(contig %in% present_in_all_countries) %>%
  left_join(., SRA_to_study, by = "SRA") %>%
  group_by(contig, study) %>%
  mutate(present = ifelse(sum(reads >0), TRUE, FALSE)) %>%
  ungroup()

present_in_dataset <- list()
for (dataset in c("Deboutte", "Busby", "Bonilla")) {
  present_in_dataset$three_way[[dataset]] <- presence_absence %>%
    filter(study == dataset) %>%
    filter(present) %>%
    select(contig) %>%
    distinct() %>%
    unlist(use.names = FALSE)
}
present_in_dataset$four_way <- c(list(BPhage = present_in_all_countries), present_in_dataset$three_way)

core_read_presence_overlap <- list()
for (thing in names(present_in_dataset)) {
  core_read_presence_overlap[[thing]] <- ggVennDiagram(present_in_dataset[[thing]]) +
    theme(legend.position = "none")
}
core_read_presence_overlap
read_counts
  
  
##### Save files

system("mkdir -p output/R/other_studies")
for (thing in names(conitg_overlap_venn)) {
  ggsave(paste0("output/R/other_studies/conitg_overlap.", thing, ".pdf"),
         conitg_overlap_venn[[thing]], width = 6, height = 6)
}
for (thing in names(core_read_presence_overlap)) {
  ggsave(paste0("output/R/other_studies/core_read_presence_overlap.", thing,".pdf"),
         core_read_presence_overlap[[thing]], width = 6, height = 6)
}
write_csv(read_counts, "output/R/other_studies/read_counts.csv")
