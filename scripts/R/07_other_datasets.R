library(ggVennDiagram)
library(patchwork)
library(tidyverse)

present_in_all_countries <- read_lines("data/core_contigs.txt")
prevalence.Bee_pools <- read.csv("output/R/prevalence/prevalence.Bee_pools.csv") %>%
  rename(bee_pool_prevalence = prevalence_abs) %>%
  select(-prevalence_prop)
prevalence.Hives <- read.csv("output/R/prevalence/prevalence.Hives.csv") %>%
  rename(hive_prevalence = prevalence_abs) %>%
  select(-prevalence_prop)
classification <- readRDS("output/R/R_variables/classification.RDS") %>%
  tibble()


##### Contig overlap
# 
# clusters <- list()
# 
# clusters$all_BPhage <- read.delim("output/bphage_and_others_clusters.tsv", header=FALSE) %>%
#   rename(representative = V1, member = V2) %>%
#   separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")
# 
# clusters$core_BPhage <- read.delim("output/bphage_core_and_others_clusters.tsv", header=FALSE) %>%
#   rename(representative = V1, member = V2) %>%
#   separate_wider_delim(member, ",", names_sep = "_", too_few = "align_start")
# 
# conitg_overlap_venn <- list()
# for (core_or_not in names(clusters)) {
#   presence_in_datasets <- clusters[[core_or_not]] %>%
#     mutate(across(
#       starts_with("member_"),
#       ~ case_when(
#         str_detect(., "NODE") ~ "BPhage",
#         str_detect(., "Deboutte") ~ "Deboutte",
#         str_detect(., "Busby") ~ "Busby",
#         str_detect(., "Bonilla") ~ "Bonilla",
#         TRUE ~ .
#       ))) %>%
#     pivot_longer(-representative, values_drop_na = TRUE, names_to = "member", values_to = "dataset") %>%
#     select(-member) %>%
#     distinct()
#   
#   genomes_in_dataset <- list()
#   for (set in c("BPhage", "Deboutte", "Busby", "Bonilla")) {
#     genomes_in_dataset[[set]] <- presence_in_datasets %>%
#       filter(dataset == set) %>%
#       select(representative) %>%
#       unlist(use.names = FALSE)
#   }
#   
#   conitg_overlap_venn[[core_or_not]] <- ggVennDiagram(genomes_in_dataset, 
#                                                       label = "count",
#                                                       label_size = 7) +
#     theme(legend.position = "none")
# }
# 
# conitg_overlap_venn

##### Mapping

horizontal_coverage_threshold = 70
mean_depth_threshold = 1

SRA_to_study <- read.delim("data/other_datasets_SRA_accessions.tsv", header=FALSE) %>%
  rename(SRA = V1, study = V2)

stats.read.other_studies.reads <- read.delim("output/other_studies/stats.read.other_studies.reads.tsv") %>%
  select(SRA, Trimmed_pairs) %>%
  left_join(., SRA_to_study, by = "SRA") %>%
  group_by(study) %>%
  summarise(trimmed_read_pairs = sum(Trimmed_pairs),
            trimmed_reads = sum(Trimmed_pairs)*2) %>%
  mutate(study = ifelse(study == "Bonilla", "Bonilla-Rosso", study))

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

print_names <- c("Deboutte" = "Deboutte",
                 "Busby" = "Busby",
                 "Bonilla" = "Bonilla-Rosso",
                 "Sbardellati" = "Sbardellati",
                 "Feng" = "Feng")
present_in_dataset <- list()
for (dataset in unique(SRA_to_study$study)) {
  print_name <- print_names[[dataset]]
  present_in_dataset[[print_name]] <- presence_absence %>%
    filter(study == dataset) %>%
    filter(present) %>%
    select(contig) %>%
    distinct() %>%
    unlist(use.names = FALSE)
}

core_read_presence_overlap <- list()
set_names <- names(present_in_dataset)
core_read_presence_overlap$venn <- ggVennDiagram(present_in_dataset, 
                                                     label = "count",
                                                     label_size = 7,
                                                     category.names = set_names,
                                                     set_size = 5.5) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_fill_gradient(low = "#8B4513", high = "#FFC300")

vennObj <- Venn(present_in_dataset)

core_read_presence_overlap$upset <- plot_upset(vennObj,
           nintersects = 15,
           sets.bar.color  = "black",
           top.bar.color = "black",
           intersection.matrix.color = "black",
           # top.bar.show.numbers = FALSE,
           top.bar.numbers.size = 7,
           # relative_height = 5,
           # sets.bar.show.numbers = TRUE
           # relative_width = 1.5
)
core_read_presence_overlap$upset

only_intersection_bars <- tibble(intersection_set = LETTERS[1:15],
       set_count = rev(c(1,1,3,3,3,4,4,6,7,8,8,10,11,13,15))) %>%
  ggplot(aes(x = intersection_set, y = set_count)) +
  geom_col(fill = "black") +
  theme_void() +
  geom_text(
    aes(label = set_count),
    vjust = -0.5,
    size  = 9
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .175))) +
  theme(
    plot.margin = margin(5,5,5,5),
    )


# All this only because they don't let me manipulate the stuff directly... 
p <- core_read_presence_overlap$upset[[1]]
for (i in seq_along(p$layers)) {
  geom_class <- class(p$layers[[i]]$geom)[1]
  
  if (geom_class == "GeomPoint") {
    # override the point size
    p$layers[[i]]$aes_params$size <- 8
    
  } else if (geom_class %in% c("GeomLine", "GeomPath")) {
    # override the line width
    p$layers[[i]]$aes_params$linewidth <- 4
  }
}

only_points <- p +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16, face = "bold"),
    plot.margin = margin(5,5,5,5)
    )

upset_patch <- only_intersection_bars / only_points + plot_layout(heights = c(2,1))

dataset_overlap <- tibble(contig = present_in_all_countries) %>%
  mutate(Bphage = TRUE,
         Deboutte = ifelse(contig %in% present_in_dataset$Deboutte, TRUE, FALSE),
         Bonilla = ifelse(contig %in% present_in_dataset$Bonilla, TRUE, FALSE),
         Busby = ifelse(contig %in% present_in_dataset$Busby, TRUE, FALSE),
         Sbardellati = ifelse(contig %in% present_in_dataset$Sbardellati, TRUE, FALSE),
         Feng = ifelse(contig %in% present_in_dataset$Feng, TRUE, FALSE),
  ) %>%
  pivot_longer(-contig) %>%
  group_by(contig) %>%
  mutate(dataset_prevalence = sum(value)) %>%
  ungroup() %>%
  pivot_wider() %>%
  relocate(dataset_prevalence, .after = last_col()) %>%
  left_join(., prevalence.Hives, by = "contig") %>%
  left_join(., prevalence.Bee_pools, by = "contig")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()
metadata <- readRDS("output/R/R_variables/metadata.RDS")

positive_pools <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "reads") %>% # it's not really reads but it doesn't matter here
  left_join(., metadata[c("Sample_ID", "Bee_pool")], by = "Sample_ID") %>%
  rename(SRA = Bee_pool) %>%
  select(-Sample_ID) %>%
  distinct() %>%
  rbind(filtered_ab_long) %>%
  filter(contig %in% present_in_all_countries) %>%
  left_join(., SRA_to_study, by = "SRA") %>%
  mutate(study = ifelse(is.na(study), "BPhage", study),
         study = ifelse(study == "Bonilla", "Bonilla-Rosso", study)) %>%
  filter(reads > 0) %>%
  group_by(study) %>%
  summarize(pools_with_core_phages = n_distinct(SRA))
  
read_counts_and_pools <- read.delim("output/bphage_viper_output/read_stats.tsv") %>%
  tibble() %>%
  filter(!str_starts(Sample, "Blank")) %>%
  mutate(Trimmed_pairs = Trimmed_R1_plus_R2 / 2) %>%
  mutate(study = "BPhage") %>%
  select(study, Trimmed_pairs) %>%
  group_by(study) %>%
  summarise(trimmed_read_pairs = sum(Trimmed_pairs),
            trimmed_reads = sum(Trimmed_pairs) * 2) %>%
  rbind(stats.read.other_studies.reads) %>%
  mutate(pools = case_when(study == "BPhage" ~ 150,
                           study == "Bonilla-Rosso" ~ 2,
                           study == "Busby" ~ 1,
                           study == "Deboutte" ~ 102,
                           study == "Sbardellati" ~ 3,
                           study == "Feng" ~ 6),
         bees_per_pool = case_when(study == "BPhage" ~ 10,
                                   study == "Bonilla-Rosso" ~ 100,
                                   study == "Busby" ~ 75,
                                   study == "Deboutte" ~ 6,
                                   study == "Sbardellati" ~ 100,
                                   study == "Feng" ~ 100)) %>%
  left_join(., positive_pools, by = "study") %>%
  mutate(sampling_time = case_when(study == "BPhage" ~ "2020",
                                   study == "Bonilla-Rosso" ~ "2015/16",
                                   study == "Busby" ~ "2020",
                                   study == "Deboutte" ~ "2012/13",
                                   study == "Sbardellati" ~ "2023",
                                   study == "Feng" ~ "2020"), # Autumn 2020, only mentioned in "Reporting summary"
         sampling_region = case_when(study == "BPhage" ~ "Europe",
                                     study == "Bonilla-Rosso" ~ "Switzerland",
                                     study == "Busby" ~ "Texas, USA",
                                     study == "Deboutte" ~ "Belgium",
                                     study == "Sbardellati" ~ "California, USA",
                                     study == "Feng" ~ "China"),
         total_core_phages = case_when(study == "BPhage" ~ sum(dataset_overlap$Bphage),
                                     study == "Bonilla-Rosso" ~ sum(dataset_overlap$Bonilla),
                                     study == "Busby" ~ sum(dataset_overlap$Busby),
                                     study == "Deboutte" ~ sum(dataset_overlap$Deboutte),
                                     study == "Sbardellati" ~ sum(dataset_overlap$Sbardellati),
                                     study == "Feng" ~ sum(dataset_overlap$Feng)
                                 )
  ) %>%
  mutate(trimmed_reads = trimmed_read_pairs*2, .after = trimmed_read_pairs) %>%
  arrange(desc(total_core_phages))

ratio_bphage_to_all_others <- 3771661288 / (sum(read_counts_and_pools$trimmed_read_pairs)-3771661288)

##### Save files

system("mkdir -p output/R/other_studies")
# for (thing in names(conitg_overlap_venn)) {
#   ggsave(paste0("output/R/other_studies/conitg_overlap.", thing, ".pdf"),
#          conitg_overlap_venn[[thing]], width = 8, height = 6)
# }
ggsave(paste0("output/R/other_studies/core_read_presence_overlap.venn.pdf"),
       core_read_presence_overlap$venn, width = 8, height = 6)
ggsave(paste0("output/R/other_studies/core_read_presence_overlap.upset.pdf"),
       core_read_presence_overlap$upset, width = 18, height = 6)

ggsave(paste0("output/R/other_studies/core_read_presence_overlap.upset.patch.pdf"),
       upset_patch, width = 14, height = 5.5)

write_csv(read_counts_and_pools, "output/R/other_studies/read_counts_and_pools.csv")
write_csv(dataset_overlap, "output/R/other_studies/dataset_overlap.csv")

# For convenience, to avoid backtracking
# new_classification_df <- dataset_overlap %>%
#   select(-Bphage) %>%
#   mutate(Prevalence_other_datasets = dataset_prevalence-1, .after = "contig") %>%
#   select(-c(dataset_prevalence, hive_prevalence, bee_pool_prevalence)) %>%
#   rename(Present_in_Deboutte = Deboutte,
#          Present_in_Bonilla = Bonilla,
#          Present_in_Busby = Busby,
#          Present_in_Sbardellati = Sbardellati,
#          Present_in_Feng = Feng) %>%
#   left_join(classification, ., by = "contig")
# write_csv(new_classification_df, "output/R/classification.csv")
# saveRDS(new_classification_df, "output/R/R_variables/classification.RDS")



