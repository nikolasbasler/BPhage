start_time <- Sys.time()
library(phyloseq)
library(furrr)
library(patchwork)
library(vegan)
library(ape)
library(ggpubr)
library(gridExtra)
library(VennDiagram)
library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(UpSetR)

source("scripts/R/helpers/tpm_and_reads_per_kb.R")
source("scripts/R/helpers/read_count_stats.R")
source("scripts/R/helpers/relative_abundance.R")
source("scripts/R/helpers/alpha_diversity.R")
source("scripts/R/helpers/beta_diversity.R")
source("scripts/R/helpers/decontam.R")
source("scripts/R/helpers/venn.R")
source("scripts/R/helpers/vcontact.R")

set.seed(1)
iterations <- 1000
# color_vector <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set1")[c(1,2,8)])
# color_vector <- colorRampPalette(color_vector)(43)
color_vector <- colorRampPalette(brewer.pal(n = 8, name = "Set3"))(43) %>%
  sample()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Read in data and generate abundance and TPM tables ####

metadata <- read.csv("data/metadata.csv") %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  mutate(Season = factor(Season, levels = c("spr", "sum", "aut"))) %>%
  mutate(Gut_part = factor(Gut_part, levels = c("mid", "ile", "rec"))) %>%
  mutate(Health = factor(Health, levels = c("very_bad", "bad", "good", "very_good")))
sample_order <- metadata %>% 
  arrange(Country, Season, Hive_ID, Gut_part) %>%
  select(Sample_ID) %>%
  unlist(use.names = FALSE)
metadata <- metadata %>%
  mutate(Sample_ID = factor(Sample_ID, levels=sample_order))
row.names(metadata) <- metadata$Sample_ID

# The inputs are exported tables from Cytoscape.
microvirus_contigs <-  list()
microvirus_contigs[["0"]] <- read_lines("data/cytoscape_yFiles_Organic_layout_bin_0.topleft_blop")
for (subgraph in as.character(1:2)) { # <- mind the as.character()
  microvirus_contigs[[subgraph]] <- read.csv(paste0("output/vcontact3/bphage_vcontact3_b38_with_inphared/graph.bin_", subgraph, ".cyjs default node.csv")) %>%
    select(X_nx_name) %>%
    unlist(use.names = FALSE)
}

classification_gnmd <- read.csv("output/R/phage.filt.gnmd.classification.csv") %>%
  mutate(contig_length = contig_length/1000) %>%
  rename(length_kb = contig_length)

present_in_all_countries <- read_lines("data/core_contigs.txt")

# classification <- read.csv("output/vcontact3/previous_runs/vcontact3/final_assignments.csv") %>%
# classification <- read.csv("output/vcontact3/previous_runs/vcontact3_with_inphared/final_assignments.csv") %>%

classification <- read.csv("output/vcontact3/bphage_vcontact3_b38_with_inphared/final_assignments.csv")  %>% 
  # filter(str_detect(GenomeName, "NODE") | str_detect(GenomeName, "Busby") | str_detect(GenomeName, "Bonilla") | str_detect(GenomeName, "Deboutte"))
  filter(str_detect(GenomeName, "NODE"))
classification <- pick_ambiguous_taxa(vcontact_output = classification,
                                      taxlevel_to_pick = "Subfamily")
classification <- pick_ambiguous_taxa(vcontact_output = classification,
                                      taxlevel_to_pick = "Genus")

classification <- classification %>%
mutate(Kingdom = "",  # No Kingdom column in vcontact's output??
         Species = "" ) %>%
  rename(contig = GenomeName,
         length_kb = Size..Kb.,
         Realm = realm..prediction.,
         Phylum = phylum..prediction.,
         Class = class..prediction.,
         Order = order..prediction.,
         Family = family..prediction.) %>%
  select(contig, length_kb, Realm, Kingdom, Phylum, Class, Order, Family, Subfamily, Genus, Species) %>%
  mutate(Family = ifelse(contig %in% microvirus_contigs[["0"]], "Microvirus_cluster_A", Family)) %>%
  mutate(Family = ifelse(contig %in% microvirus_contigs[["1"]], "Microvirus_cluster_B", Family)) %>%
  mutate(Family = ifelse(contig %in% microvirus_contigs[["2"]], "Microvirus_cluster_C", Family)) %>%
  mutate(Order = ifelse(str_detect(Family, "Microvirus"), "Microviruses", Order)) %>%
  mutate(lowest_taxon = apply(., 1, function(row) tail(row[row != ""], 1))) %>%
  mutate(across(everything(), ~ifelse(. == "", "Unclassified", .))) %>%
  inner_join(., classification_gnmd[c("contig", "provirus", "proviral_length", 
                                     "gene_count", "viral_genes", "host_genes",
                                     "checkv_quality", "miuvig_quality", 
                                     "completeness", "completeness_method",
                                     "contamination", "kmer_freq", "warnings")],
            by = "contig") %>%
  mutate(predicted_genome_length = length_kb / (completeness/100)) %>%
  mutate(Host_group = "all") %>% 
  mutate(Core = ifelse(contig %in% present_in_all_countries, "yes", "no"))

# This would change each "Unclassified" entry to unclassified_<taxlevel>_of_<higher_taxlevel>
# Disabled because it would break things downstream.
# taxlevels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species")
# higher_tax <- list("Kingdom" = "Realm", 
#                    "Phylum" = "Kingdom", 
#                    "Class" = "Phylum",
#                    "Order"= "Class",
#                    "Family" = "Order",
#                    "Subfamily" = "Family",
#                    "Genus" = "Subfamily",
#                    "Species" = "Genus")
# for (tl in taxlevels) {
#   classification <- classification %>%
#     mutate(!!tl := ifelse(.data[[tl]] == "Unclassified", paste0("unclassified_", tl, "_of_", .data[[higher_tax[[tl]]]]), .data[[tl]]))
#   
# }

contig_order <- classification %>%
  arrange(Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  select(contig) %>%
  unlist(use.names = FALSE)

classification <- classification %>%
  mutate(contig = factor(contig, levels=contig_order))

phage_abundance <- read.csv("output/R/phage.filt.abundance.contig.csv") %>%
  select(!contains("Blank"))

# ANI
bphage_and_others_ani.tsv <- read.delim("output/ani/bphage_and_others_ani.tsv.gz")

classification %>% 
  mutate(classified_on_family_level = ifelse(Family == "Unclassified", "no", "yes")) %>%
  group_by(classified_on_family_level) %>%
  count()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Taxon and core pies ####

core_pie <- list()
top_dozen_taxes <- list()
taxlevels <- c("Order", "Family")
for (tl in taxlevels) {
  pie_tibble <- classification %>%
    filter(contig %in% present_in_all_countries) %>%
    select(all_of(tl)) %>%
    table() %>%
    as.data.frame() %>%
    # rename(Lowest_taxon = Var1, Genomes = Freq) %>%
    arrange(desc(Freq)) %>%
    mutate(Tax = factor(.data[[tl]], levels = .data[[tl]]))
  
  top_dozen_taxes[[tl]] <- pie_tibble %>%
    head(n=12) %>%
    select(all_of(tl)) %>%
    unlist(use.names = FALSE) %>%
    factor()
  
  core_pie[[tl]] <- pie_tibble %>%
    ggplot(aes(x = "", y=Freq, fill = Tax)) +
    geom_bar(stat = "identity", color="black", linewidth=0.2) +
    coord_polar("y") +
    theme_classic() +
    labs(x = NULL, y = NULL) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_brewer(palette="Set3", name = tl)
}

taxlevels <- c("Order", "Family")
taxon_pie <- list()
taxon_pie_overlap <- list()
for (tl in taxlevels) {
  set3_colors <- brewer.pal(n = 12, name = "Set3")
  custom_colors <- setNames(set3_colors, top_dozen_taxes[[tl]])
  
  pie_tible <- table(classification[[tl]]) %>%
    as.data.frame() %>%
    rename(Taxon = Var1, Genomes = Freq) %>%
    mutate(Taxon = as.character(Taxon),
           Taxon = ifelse(Genomes < 10, "Other (<10)", Taxon)) %>%
    group_by(Taxon) %>%
    mutate(Genomes = sum(Genomes)) %>%
    distinct() %>%
    ungroup() %>%
    arrange(desc(Genomes)) %>%
    mutate(Taxon = factor(Taxon, levels = Taxon,))
  
  taxon_pie[[tl]] <- pie_tible %>%
    ggplot(aes(x = "", y=Genomes, fill = Taxon)) +
    geom_bar(stat = "identity", color="black", linewidth=0.2) +
    coord_polar("y") +
    theme_classic() +
    labs(x = NULL, y = NULL) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_brewer(palette="Set3")
  
  taxon_pie_overlap[[tl]] <- taxon_pie[[tl]] +
    scale_fill_manual(values = custom_colors, na.value = "white")
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Generate abundance and TPM tables ####

## Generate abundance tables for different taxonomic levels.####
taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm", "Host_group", "Core")
phage_ab <- list()
phage_ab <- taxlevels %>%
  set_names() %>%
  map(~tax_sum(., ab_table=phage_abundance,
               classif = classification))
phage_lengths <- list()
phage_lengths <- taxlevels %>%
  set_names() %>%
  map(~tax_lengths(., classif = classification))

taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
hostgroups <- classification$Host_group %>% unique()
phage_ab_hostgroup <- list()
for (hostg in hostgroups) {
  hostg_filt <- hostg_filter(hg = hostg, ab_table = phage_ab$contig,
          classif = classification)
  phage_ab_hostgroup[[hostg]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=hostg_filt,
                 classif = classification))
}

taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
phage_ab_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  core_filt <- core_filter(core_y_n = core_or_not, 
                             ab_table = phage_ab$contig,
                             classif = classification)
  phage_ab_core_or_not[[core_or_not]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=core_filt,
                 classif = classification))
}

## Generate contig abundance table with merged metadata variables 
# (e.g. all gut parts or all samples from the same country merged). Used 
# for prevalence plots
meta_merges <- list(Bee_pools = "Gut_part", 
                    Countries = c("Season", "Gut_part", "Hive_ID"),
                    Hives = c("Season", "Gut_part"))
phage_ab_meta_merges <- list()
for (merge in names(meta_merges)) {
  phage_ab_meta_merges[[merge]] <- meta_sum(ab_table = phage_ab$contig, 
                                            meta_vars = meta_merges[[merge]])
}

## Generate TPM tables for different taxonomic levels. ####
phage_tpm <- list()
for (lvl in names(phage_ab)) {
  phage_tpm[[lvl]] <- calc_tpm(abtable = phage_ab[[lvl]], 
                               level = lvl, 
                               lengths_df = phage_lengths[[lvl]])
}

phage_tpm_hostgroup <- list()
for (hostg in hostgroups) {
  for (lvl in names(phage_ab_hostgroup[[hostg]])) {
    phage_tpm_hostgroup[[hostg]][[lvl]] <- calc_tpm(abtable = phage_ab_hostgroup[[hostg]][[lvl]], 
                                                    level = lvl, 
                                                    lengths_df = phage_lengths[[lvl]])
  }
}

phage_tpm_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  for (lvl in names(phage_ab_core_or_not[[core_or_not]])) {
    phage_tpm_core_or_not[[core_or_not]][[lvl]] <- calc_tpm(abtable = phage_ab_core_or_not[[core_or_not]][[lvl]], 
                                                            level = lvl, 
                                                            lengths_df = phage_lengths[[lvl]])
  }
}

## Generate absolute viral load tables for different taxonomic levels
samples_with_vlp_counts <- metadata %>%
  filter(!is.na(VLPs_per_ul)) %>%
  rownames()
viral_loads <- phage_tpm$contig %>%
  select(contig, any_of(samples_with_vlp_counts)) %>% 
  filter(!if_all(-contig, ~. == 0)) %>% 
  pivot_longer(-contig, names_to = "Sample_ID") %>%
  left_join(., metadata[c("Sample_ID", "VLPs_per_ul")], by="Sample_ID") %>%
  mutate(viral_load = value * VLPs_per_ul) %>%
  select(-c(value, VLPs_per_ul)) %>%
  pivot_wider(names_from = Sample_ID, values_from = viral_load) 

phage_load <- list()
for (lvl in names(phage_ab)) {
  phage_load[[lvl]] <- tax_sum(ab_table = viral_loads, tax_level = lvl, classif = classification)
}

taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
phage_load_hostgroup <- list()
for (hostg in hostgroups) {
  hostg_filt <- hostg_filter(hg = hostg, ab_table = phage_load$contig,
                             classif = classification)
  phage_load_hostgroup[[hostg]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=hostg_filt,
                 classif = classification))
}

taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
phage_load_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  core_filt <- core_filter(core_y_n = core_or_not, 
                           ab_table = phage_load$contig,
                           classif = classification)
  phage_load_core_or_not[[core_or_not]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=core_filt,
                 classif = classification))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# ANI comparison ####

ANI_bphage <- bphage_and_others_ani.tsv %>%
  # select(id1, id2, tani) %>%
  # filter(id1 > id2) %>%
  filter(str_detect(id1, "NODE") & str_detect(id2, "NODE"))

aANI_bphage <- ANI_bphage %>%
  select(tani) %>%
  mutate(subset = "all")
aANI_core <- ANI_bphage %>%
  filter(id1 %in% present_in_all_countries & id2 %in% present_in_all_countries) %>%
  select(tani) %>%
  mutate(subset = "core")
aANI_non_core <- ANI_bphage %>%
  filter((!id1 %in% present_in_all_countries) & (!id2 %in% present_in_all_countries)) %>%
  select(tani) %>%
  mutate(subset = "non-core")

aANI_df <- rbind(aANI_bphage, aANI_core, aANI_non_core)
kruskal_results <- aANI_df %>%
  summarize(pvalue = kruskal.test(tani~subset)$p.value,
            test_stat = kruskal.test(tani~subset)$statistic,
            deg_freedom = kruskal.test(tani~subset)$parameter)
aANI_stats <- aANI_df %>%
  group_by(subset) %>%
  summarise(mean_ani = mean(tani),
            median_ani = median(tani),
            IQR_ani = IQR(tani),
            sd_ani = sd(tani),
            min_ani = min(tani),
            max_ani = max(tani))

aANI_boxplot <- aANI_df %>%
  # sample_frac(0.001) %>% # For fast testing.
  ggplot(aes(x = subset, y = tani)) +
  geom_boxplot() +
  # geom_boxplot(outlier.shape = NA) +
  geom_pwc(method="wilcox.test", label="p.adj.signif",
           p.adjust.method="BH", hide.ns = TRUE) + #, bracket.nudge.y = -0.8) + #,
           # y.position = c(0.1, 0.13, 0.16)) +
  # scale_y_continuous(limits = c(0, 0.25)) +
  # scale_y_continuous(limits = c(0, quantile(aANI_df$tani, 0.99))) +
  labs(title = paste0("Krusil-Wallis H = ", round(kruskal_results$test_stat), ", p = ",kruskal_results$pvalue))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Completeness and genome length comparison core/non-core####

completeness_and_genome_length <- list()
completeness_and_genome_length$completeness <- classification %>%
  mutate(completeness = ifelse(contig == "NODE_A185_length_1971_cov_16.487328_RO_26022_aut_mid_d", NA, completeness)) %>% # That's the one Picobirna contig that got a completeness estimate from CheckV.
  ggplot(aes(x = Core, y = completeness)) +
  geom_boxplot() +
  geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE) +
  labs(title = "Completeness")
completeness_and_genome_length$predicted_length <- classification %>%
  ggplot(aes(x = Core, y = predicted_genome_length)) +
  geom_boxplot() +
  geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE) +
  labs(title = "Predicted genome length")
completeness_and_genome_length$contig_length <- classification %>%
  ggplot(aes(x = Core, y = length_kb)) +
  geom_boxplot() +
  geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE) +
  labs(title = "Assembled contig length")

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Set min seq thresholds for rarefaction. ####

# 50% horizontal coverage cutoff
# count_stats <- report_stats(df = phage_abundance,
#                             thresholds=c(2689, 1002)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
# min_seq_count <- 2713 # For quantile_90_ratio < 1000. 11 discarded samples (10 mid, 1 ile).

# 70% horizontal coverage cutoff
count_stats <- report_stats(df = phage_abundance,
                            thresholds=c(2437, 775)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
count_stats
# min_seq_count <- 2442 # For quantile_90_ratio < 1000. 12 discarded samples (9 mid, 3 ile).
min_seq_count <- 1548 # Looks like a natural breaking point. 4 discarded samples (all mid).
discarded <- discards(count_stats$ratios, min_seq_count)$discarded
lost_bees <- discards(count_stats$ratios, min_seq_count)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

min_seq_count_core_or_not <- list()
count_stats_core_or_not <- list()
count_stats_core_or_not$no <- report_stats(df = phage_ab_core_or_not$no$contig,
                            thresholds=c(961, 7356))
count_stats_core_or_not$no
min_seq_count_core_or_not$no <- 961 # Looks like a natural breaking point. 7 discarded samples (5 mid, 2 ile).
discarded <- discards(count_stats_core_or_not$no$ratios, min_seq_count_core_or_not$no)$discarded
lost_bees <- discards(count_stats_core_or_not$no$ratios, min_seq_count_core_or_not$no)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

count_stats_core_or_not$yes <- report_stats(df = phage_ab_core_or_not$yes$contig,
                                     thresholds=c(462, 2375))
count_stats_core_or_not$yes
min_seq_count_core_or_not$yes <- 870 # Looks like a natural breaking point. 1 discarded ile sample.
discarded <- discards(count_stats_core_or_not$yes$ratios, min_seq_count_core_or_not$yes)$discarded
lost_bees <- discards(count_stats_core_or_not$yes$ratios, min_seq_count_core_or_not$yes)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Alpha ####
# Full set and core or not
alpha_start <- Sys.time()
met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha <- list()
alpha_core_or_not <- list()
for (tlvl in taxlevels) {
  alpha[[tlvl]] <- alpha_stats(df = phage_ab[[tlvl]], 
                               meta_vars = met_v, 
                               min_seq = min_seq_count,
                               df_lengths = phage_lengths[[tlvl]])
  for (core_or_not in unique(classification$Core)) {
    alpha_core_or_not[[core_or_not]][[tlvl]] <- alpha_stats(df = phage_ab_core_or_not[[core_or_not]][[tlvl]], 
                                 meta_vars = met_v, 
                                 min_seq = min_seq_count_core_or_not[[core_or_not]],
                                 df_lengths = phage_lengths[[tlvl]])
  }
}

# By country and core or not by country
met_v <- c("Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_by_country <- list()
alpha_by_country_core_or_not <- list()
for (countr in levels(metadata$Country)) {
  for (tlvl in taxlevels) {
    count_filt_ab <- phage_ab[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    alpha_by_country[[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                 meta_vars = met_v, 
                                 min_seq = min_seq_count,
                                 df_lengths = phage_lengths[[tlvl]])
    for (core_or_not in unique(classification$Core)) {
      count_filt_ab <- phage_ab_core_or_not[[core_or_not]][[tlvl]] %>%
        select(all_of(tlvl), starts_with(countr))
      alpha_by_country_core_or_not[[core_or_not]][[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                              meta_vars = met_v, 
                                                              min_seq = min_seq_count_core_or_not[[core_or_not]],
                                                              df_lengths = phage_lengths[[tlvl]])
    }
  }
}

# Absolute counts full (measured) set and core or not
met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_abs <- list()
alpha_abs_core_or_not <- list()
for (tlvl in taxlevels) {
  alpha_abs[[tlvl]] <- alpha_stats(df = phage_load[[tlvl]], 
                                   absolut_values = TRUE,
                                   meta_vars = met_v)
  for (core_or_not in unique(classification$Core)) {
    alpha_abs_core_or_not[[core_or_not]][[tlvl]] <- alpha_stats(df = phage_load_core_or_not[[core_or_not]][[tlvl]], 
                                     absolut_values = TRUE,
                                     meta_vars = met_v)
  }
}

# Absolute counts by country and core or not by country
met_v <- c("Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
alpha_abs_by_country <- list()
alpha_abs_by_country_core_or_not <- list()
for (countr in levels(metadata$Country)) {
  for (tlvl in taxlevels) {
    count_filt_ab <- phage_load[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    alpha_abs_by_country[[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                     absolut_values = TRUE,
                                     meta_vars = met_v)
    for (core_or_not in unique(classification$Core)) {
      count_filt_ab <- phage_load_core_or_not[[core_or_not]][[tlvl]] %>%
        select(all_of(tlvl), starts_with(countr))
      alpha_abs_by_country_core_or_not[[core_or_not]][[countr]][[tlvl]] <- alpha_stats(df = count_filt_ab, 
                                                                  absolut_values = TRUE,
                                                                  meta_vars = met_v)
    }
  }
}
alpha_end <- Sys.time()
alpha_end - alpha_start
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Beta ####
# Full set
beta_start <- Sys.time()
met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_dist <- list()
beta_plot_list <- list()
for (tlvl in taxlevels) {
  beta_dist[[tlvl]] <- ordination(df = phage_ab[[tlvl]], 
                                        meta_vars = met_v, 
                                        min_seq = min_seq_count,
                                        df_lengths = phage_lengths[[tlvl]])
  
  beta_plot_list[[tlvl]] <- beta_plot(ordination_list = beta_dist[[tlvl]]$ord_list, 
                                      meta_vars = met_v,
                                      mapped_reads = count_stats$ratios)
}

# Core or not
met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_dist_core_or_not <- list()
beta_plot_list_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  for (tlvl in taxlevels) {
    beta_dist_core_or_not[[core_or_not]][[tlvl]] <- ordination(df = phage_ab_core_or_not[[core_or_not]][[tlvl]], 
                                    meta_vars = met_v, 
                                    min_seq = min_seq_count,
                                    df_lengths = phage_lengths[[tlvl]])
    
    beta_plot_list_core_or_not[[core_or_not]][[tlvl]] <- beta_plot(ordination_list = beta_dist_core_or_not[[core_or_not]][[tlvl]]$ord_list, 
                                        meta_vars = met_v,
                                        mapped_reads = count_stats_core_or_not[[core_or_not]]$ratios)
  }
}

# Absolute counts
met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_abs_dist <- list()
beta_abs_plot_list <- list()
for (tlvl in taxlevels) {
  beta_abs_dist[[tlvl]] <- ordination(df = phage_load[[tlvl]],
                                      meta_vars = met_v,
                                      absolute_values = TRUE)
  beta_abs_plot_list[[tlvl]] <- beta_plot(beta_abs_dist[[tlvl]]$ord_list,
                                          meta_vars = met_v,
                                          mapped_reads = count_stats$ratios)

}

# Absolute counts core or not
met_v <- c("Country", "Season", "Health")
taxlevels <- c("contig", "Genus", "Family")
beta_abs_dist_core_or_not <- list()
beta_abs_plot_list_core_or_not <- list()
for (core_or_not in unique(classification$Core)) {
  for (tlvl in taxlevels) {
    beta_abs_dist_core_or_not[[core_or_not]][[tlvl]] <- ordination(df = phage_load_core_or_not[[core_or_not]][[tlvl]],
                                        meta_vars = met_v,
                                        absolute_values = TRUE)
    beta_abs_plot_list_core_or_not[[core_or_not]][[tlvl]] <- beta_plot(beta_abs_dist[[tlvl]]$ord_list,
                                            meta_vars = met_v,
                                            mapped_reads = count_stats_core_or_not[[core_or_not]]$ratios)
    
  }
}
beta_end <- Sys.time()
beta_end - beta_start
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Heatmaps ####

order_heatmaps_row <- list()
plotted_samples <- list()
plotted_contigs <- list()
for (ord in unique(classification$Order)) {
  if (ord == "Unclassified") {
    next
  }
    ord_tpm <- phage_tpm$contig %>%
      inner_join(., classification, by="contig") %>%
      filter(Order==ord) %>%
      select(colnames(phage_tpm$contig)) %>%
      select(which(colSums(. != 0) > 0))
    plotted_samples[[ord]] <- ncol(ord_tpm)
    plotted_contigs[[ord]] <- nrow(ord_tpm)
    order_heatmaps_row[[ord]] <- contig_heatmap(df = ord_tpm,
                                                   classif = classification)
}

tax_collapse_heatmaps <- list()
plotted_samples_tax_collapse <- list()
plotted_contigs_tax_collapse <- list()
taxlevels <- c("Class", "Order", "Family")
lower_tax <- list("Class" = "Order", "Order" = "Family", "Family" = "Genus")
for (tl in taxlevels) {
  lower_tl <- lower_tax[[tl]]
  for (tax in unique(classification[[tl]])) {
    lower_taxes <- classification %>%
      filter(.data[[tl]] == tax) %>%
      select(all_of(lower_tl)) %>%
      distinct() %>%
      unlist(use.names = FALSE)
    
    plot_tbl <- phage_tpm[[lower_tl]] %>%
    pivot_longer(-all_of(lower_tl), names_to = "sample", values_to = "TPM") %>% 
    filter(.data[[lower_tl]] %in% lower_taxes) %>%
    filter(TPM >0) %>%
    inner_join(., metadata, by = join_by(sample == Sample_ID))
    
    plotted_samples_tax_collapse[[tl]][[tax]] <- unique(plot_tbl$sample) %>% length()
    
    plotted_contigs_tax_collapse[[tl]][[tax]] <- unique(plot_tbl[[lower_tl]]) %>% length()
    
    tax_collapse_heatmaps[[tl]][[tax]] <- plot_tbl %>%
    ggplot(aes(x = sample, y = .data[[lower_tl]], fill = TPM)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_row(~Country+Season, scales="free_x", space="free") +
    scale_fill_gradient(limits = c(0, 1)) +
    labs(title = paste0(tl, " - ", tax))
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Host group histrograms: ####
hostgroup_hist <- list()
for (hostgroup in phage_tpm$Host_group$Host_group) {
  hostgroup_hist[[hostgroup]] <- phage_tpm$Host_group %>%
    pivot_longer(-Host_group) %>%
    filter(Host_group==hostgroup) %>%
    ggplot(aes(x=value)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("TPMs of hostgroup ",hostgroup))
}

# Core or not histograms: ####
core_or_not_hist <- list()
for (core_or_not in phage_tpm$Core$Core) {
  core_or_not_hist[[core_or_not]] <- phage_tpm$Core %>%
    pivot_longer(-Core) %>%
    filter(Core==core_or_not) %>%
    ggplot(aes(x=value)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("TPMs of core \"",core_or_not, "\""))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Average TPM ####
met_v <- c("Country", "Season", "Gut_part", "Health", "Sample_ID")
tax_levels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
average_tpm <- list()
average_tpm$Host_group <- average_tpm_bar_plot(tpm_table = phage_tpm$Host_group,
                                               tl="Host_group", 
                                               hg="Host_group",
                                               meta_vars=met_v,
                                               hg_or_core = "Host_group")
for (hgr in names(phage_tpm_hostgroup)) {
  for (tlvl in tax_levels) {
    average_tpm[[hgr]][[tlvl]] <- average_tpm_bar_plot(
      tpm_table = phage_tpm_hostgroup[[hgr]][[tlvl]],
      tl = tlvl,
      hg = hgr,
      meta_vars = met_v,
      threshold_for_other=0.01,
      hg_or_core = "Host_group")
  }
}

average_tpm_core_or_not <- list()
average_tpm_core_or_not$Core_or_not <- average_tpm_bar_plot(tpm_table = phage_tpm$Core,
                                                            tl = "Core_or_not",
                                                            hg = "Core_or_not",
                                                            meta_vars = met_v,
                                                            hg_or_core = "Core?",
                                                            threshold_for_other = 0)
for (core_or_not in names(phage_tpm_core_or_not)) {
  for (tlvl in tax_levels) {
    average_tpm_core_or_not[[core_or_not]][[tlvl]] <- average_tpm_bar_plot(
      tpm_table = phage_tpm_core_or_not[[core_or_not]][[tlvl]],
      tl = tlvl,
      hg = core_or_not,
      meta_vars = met_v,
      threshold_for_other=0.01,
      hg_or_core = "Core?")
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Core TPM ####

metavars <- c("Season", "Gut_part", "Health", "Country")
core_tpm_stats <- list()
core_tpm_plots <- list()
for (mvar in metavars) {
  core_tpm_tbl <- phage_tpm$contig %>%
    filter(contig %in% present_in_all_countries) %>%
    pivot_longer(-contig, names_to = "Sample_ID", values_to = "core_TPM") %>%
    left_join(., metadata, by = "Sample_ID") %>%
    select(all_of(mvar), core_TPM) %>%
    filter(core_TPM > 0 )
  
  kruskal_results <- core_tpm_tbl %>%
    summarize(pvalue = kruskal.test(core_TPM~.data[[mvar]])$p.value,
              test_stat = kruskal.test(core_TPM~.data[[mvar]])$statistic,
              deg_freedom = kruskal.test(core_TPM~.data[[mvar]])$parameter)
  
  core_tpm_stats[[mvar]] <- core_tpm_tbl %>%
    group_by(.data[[mvar]]) %>%
    summarise(mean = mean(core_TPM),
              median = median(core_TPM),
              min = min(core_TPM),
              max = max(core_TPM),
              sd = sd(core_TPM),
              IQR = IQR(core_TPM))
  
  core_tpm_plots[[mvar]] <- core_tpm_tbl %>%
    ggplot(aes(x = .data[[mvar]], y = core_TPM)) +
    geom_boxplot() +
    geom_pwc(method="wilcox.test", label="p.adj.signif",
             p.adjust.method="BH", hide.ns = TRUE) +
    labs(title = paste0(mvar, " - KW: H = ", round(kruskal_results$test_stat), " p = ", round(kruskal_results$pvalue)))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalence ####
prevalence_histo <- list()
for (merge in names(phage_ab_meta_merges)) {
  prevalence_histo[[merge]] <- prevalence_histogram(abtable = phage_ab_meta_merges[[merge]],
                                                          plot_title = merge)
}

# present_in_all_countries <- prevalence_histo$Countries$table %>%
#   filter(prevalence_prop ==1) %>%
#   select(group) %>%
#   unlist(use.names = FALSE) # %>%
#   # write_lines("output/temp.txt")


# 
# met_v <- c("Country", "Season", "Gut_part", "Health")
# tax_levels <- c("contig", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
# prevalence_plots <- list()
# for (hgr in names(phage_tpm_hostgroup)) {
#   for (tlvl in tax_levels) {
#     # trsh <- 3
#     if (tlvl=="contig") {
#       trsh <- 20
#     } else {
#       trsh <- 1
#     }
#     prevalence_plots[[hgr]][[tlvl]] <- prevalence_bar_plot(
#       abtable = phage_tpm_hostgroup[[hgr]][[tlvl]],
#       tl = tlvl, 
#       hg = hgr,
#       meta_vars = met_v,
#       threshold_for_other = trsh)
#   }
# }



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Country by country average TPM and prevalence: ####
tax_levels <- c("Family", "Order")
met_v <- c("Season", "Gut_part", "Health", "Sample_ID")
average_tpm_per_country <- list()
average_tpm_per_country_metavar <- list()
prevalence_plots_per_country <- list()
for (countr in levels(metadata$Country)) {
  tpm_filt <- phage_tpm$Core %>%
    select("Core", starts_with(countr))
  average_tpm_per_country_metavar[[countr]] <- average_tpm_bar_plot(tpm_table = tpm_filt,
                                                                        tl = "Core_or_not",
                                                                        hg = "Core_or_not",
                                                                        meta_vars = met_v,
                                                                        hg_or_core = "Core",
                                                                        threshold_for_other = 0)
  for (tlvl in tax_levels) {
    for (core_or_not in names(phage_tpm_core_or_not)) {
      tpm_filt <- phage_tpm_core_or_not[[core_or_not]][[tlvl]] %>%
          select(all_of(tlvl), starts_with(countr))
      country_tpm_plots <- average_tpm_bar_plot(
        tpm_table = tpm_filt,
        tl = tlvl, 
        hg = core_or_not, 
        meta_vars = met_v,
        title_prefix = paste0(countr, " - "),
        threshold_for_other = 0.03, # <==== MIND THIS!
        hg_or_core = "Core")$plots 
      average_tpm_per_country[[core_or_not]][[tlvl]][[countr]] <- 
          country_tpm_plots$Sample_ID /
          (country_tpm_plots$Season + country_tpm_plots$Gut_part + country_tpm_plots$Health)
        
        ab_filt <- phage_ab_core_or_not[[core_or_not]][[tlvl]] %>%
            select(all_of(tlvl), starts_with(countr))
        number_of_samples <- ncol(ab_filt)-1
        prevalence_plots_per_country[[core_or_not]][[tlvl]][[countr]] <- prevalence_bar_plot(
          abtable = ab_filt,
          tl = tlvl,
          hg = core_or_not,
          meta_vars = c("Season", "Gut_part", "Health"),
          title_prefix = paste0(countr, " (",number_of_samples," samples) - "),
          threshold_for_other=3)
    }
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Upset plot for country-wide sharing ####
# A Venn diagram wouldn't be feasible here because there are 8 countries. The 
# Plot will be saved into the venn folder though
country_upset_plot <- upset_country(abtable = phage_ab_meta_merges$Countries)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Genera per family in core vs non-core ####

genera_per_family_tbl <- classification %>%
  group_by(Core, Family) %>%
  summarise(genera_per_famliy = n_distinct(Genus),
            .groups = "drop")

genera_per_family_plot <- genera_per_family_tbl %>%
  ggplot(aes(x = Core, y = genera_per_famliy)) +
  geom_boxplot() +
  # geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE)
  geom_pwc(method="wilcox.test", hide.ns = TRUE)

genera_per_family_stats <- genera_per_family_tbl %>%
  group_by(Core) %>%
  summarise(mean = mean(genera_per_famliy),
            median = median(genera_per_famliy),
            sd = sd(genera_per_famliy), 
            max = max(genera_per_famliy),
            min = min(median(genera_per_famliy)),
            iqr = IQR(genera_per_famliy))

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalenve Venn diagrams ####
# Output written to files already here.
venn_stats <- list()
tax_levels <- c("contig", "Genus", "Family")
met_v <- c("Season", "Gut_part", "Health")
for (tlvl in tax_levels) {
  for (m_var in met_v) {
    system(paste0("mkdir -p output/R/venns/core_yes/", tlvl, "/", m_var))
    system(paste0("mkdir -p output/R/venns/core_no/", tlvl, "/", m_var))
    system(paste0("mkdir -p output/R/venns/all/", tlvl, "/", m_var))
    for (cntr in levels(metadata$Country)) {
      venn_stats$all[[tlvl]][[m_var]][[cntr]] <- prevalence_venn(
        abtable = phage_ab_core_or_not$yes[[tlvl]],
        meta_var = m_var,
        country = cntr,
        subset = "all")
      for (core_or_not in names(phage_tpm_core_or_not)) {
        venn_stats[[core_or_not]][[tlvl]][[m_var]][[cntr]] <- prevalence_venn(
          abtable = phage_ab_core_or_not$yes[[tlvl]],
          meta_var = m_var,
          country = cntr,
          subset = paste0("core_",core_or_not))
        }
    }
  }
}

here_time <- Sys.time()
here_time - start_time
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Save files ####

## Metadata and classification variables
system("mkdir -p output/R/R_variables")
saveRDS(metadata, "output/R/R_variables/metadata.RDS")
saveRDS(classification, "output/R/R_variables/classification.RDS")
saveRDS(classification_gnmd, "output/R/R_variables/classification_gnmd.RDS")

## Taxon and core pie
system("mkdir -p output/R/taxon_pies")
for (taxlevel in names(core_pie)) {
  ggsave(paste0("output/R/taxon_pies/core_pie.", taxlevel, ".pdf"),
         core_pie[[taxlevel]], width = 12, height = 10)
}
for (tax in names(taxon_pie)) {
  ggsave(paste0("output/R/taxon_pies/taxon_pie.", tax,".pdf"), taxon_pie[[tax]], width = 12, height = 10)
}
for (tax in names(taxon_pie_overlap)) {
  ggsave(paste0("output/R/taxon_pies/taxon_pie_overlap.", tax,".pdf"), taxon_pie_overlap[[tax]], width = 8, height = 6)
}

## aANI boxplot
ggsave("output/R/aANI_boxplot.pdf", aANI_boxplot, width = 6, height = 6)
write_csv(aANI_stats, "output/R/aANI_stats.csv")

## Completeness and genome length boxplots
system("mkdir -p output/R/completeness_and_genome_length")
for (metric in names(completeness_and_genome_length)) {
  ggsave(paste0("output/R/completeness_and_genome_length/", metric, ".pdf"), 
         completeness_and_genome_length[[metric]], width = 4, height = 4)
}

## Alpha diversity plots and tables ####
system("mkdir -p output/R/alpha")
system("mkdir -p output/R/alpha/alpha_all")
system("mkdir -p output/R/alpha/alpha_core_or_not")
for (tax in names(alpha)) {
  ggsave(paste0("output/R/alpha/alpha_all/alpha.",tax,".pdf"),
         alpha[[tax]]$plot, width = 12, height=10)
  write_csv(alpha[[tax]]$table,
            paste0("output/R/alpha/alpha_all/alpha.",tax,".csv"))
    for (core_or_not in names(alpha_core_or_not)) {
    ggsave(paste0("output/R/alpha/alpha_core_or_not/alpha_core.",core_or_not,".",tax,".pdf"),
           alpha_core_or_not[[core_or_not]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_core_or_not[[core_or_not]][[tax]]$table,
              paste0("output/R/alpha/alpha_core_or_not/alpha_core.",core_or_not,".",tax,".csv"))
  }
}

# By country
system("mkdir -p output/R/alpha/alpha_by_country")
system("mkdir -p output/R/alpha/alpha_by_country_core_or_not")
for (country in names(alpha_by_country)) {
  for (tax in names(alpha_by_country[[country]])) {
    ggsave(paste0("output/R/alpha/alpha_by_country/alpha_by_country.",country,".",tax,".pdf"),
           alpha_by_country[[country]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_by_country[[country]][[tax]]$table,
              paste0("output/R/alpha/alpha_by_country/alpha_by_country.",country,".",tax,".csv"))
    for (core_or_not in names(alpha_core_or_not)) {
      ggsave(paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_by_country_core.",core_or_not,".",country,".",tax,".pdf"),
             alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$plot, width = 12, height=10)
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$table,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_by_country_core.",core_or_not,".",country,".",tax,".csv"))
    }
  }
}

# Absolute counts
for (tax in names(alpha_abs)) {
  ggsave(paste0("output/R/alpha/alpha_all/alpha_abs.",tax,".pdf"),
         alpha_abs[[tax]]$plot, width = 12, height=10)
  write_csv(alpha_abs[[tax]]$table,
            paste0("output/R/alpha/alpha_all/alpha_abs.",tax,".csv"))
  for (core_or_not in names(alpha_abs_core_or_not)) {
    ggsave(paste0("output/R/alpha/alpha_core_or_not/alpha_abs_core.",core_or_not,".",tax,".pdf"),
           alpha_abs_core_or_not[[core_or_not]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_abs_core_or_not[[core_or_not]][[tax]]$table,
              paste0("output/R/alpha/alpha_core_or_not/alpha_abs_core.",core_or_not,".",tax,".csv"))
  }
}

# Absolute counts by country
for (country in names(alpha_abs_by_country)) {
  for (tax in names(alpha_abs_by_country[[country]])) {
    ggsave(paste0("output/R/alpha/alpha_by_country/alpha_abs_by_country.",country,".",tax,".pdf"),
           alpha_abs_by_country[[country]][[tax]]$plot, width = 12, height=10)
    write_csv(alpha_abs_by_country[[country]][[tax]]$table,
              paste0("output/R/alpha/alpha_by_country/alpha_abs_by_country.",country,".",tax,".csv"))
    for (core_or_not in names(alpha_abs_by_country_core_or_not)) {
      ggsave(paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_abs_by_country_core.",core_or_not,".",tax,".pdf"),
             alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$plot, width = 12, height=10)
      write_csv(alpha_by_country_core_or_not[[core_or_not]][[country]][[tax]]$table,
                paste0("output/R/alpha/alpha_by_country_core_or_not/alpha_abs_by_country_core.",core_or_not,".",tax,".csv"))
    }
  }
}

## Beta diversity plots and tables ####
system("mkdir -p output/R/beta/beta_all")
system("mkdir -p output/R/beta/beta_core_or_not")
for (tax in names(beta_plot_list)) {
  pcoa_path <- paste0("output/R/beta/beta_all/", tax, "_pcoa")
  pcoa_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_pcoa")
  hist_path <- paste0("output/R/beta/beta_all/", tax, "_hist")
  hist_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_hist")
  system(paste0("mkdir -p ", pcoa_path))
  system(paste0("mkdir -p ", hist_path))
  system(paste0("mkdir -p ", pcoa_path_core_or_not))
  system(paste0("mkdir -p ", hist_path_core_or_not))
  for (p in names(beta_plot_list[[tax]])) {
    if (p=="core") {
      width <- 17
      height <- 15
    } else {
      width <- 17
      height <- 5
    }
    for (q in names(beta_plot_list[[tax]][[p]])) {
      ggsave(paste0(pcoa_path,"/beta.", tax,".", p, ".", q, ".pcoa.pdf"),
             beta_plot_list[[tax]][[p]][[q]], width = width, height = height)
      for (core_or_not in names(beta_plot_list_core_or_not)) {
        ggsave(paste0(pcoa_path_core_or_not,"/beta_core.", core_or_not, ".", tax,".", p, ".", q, ".pcoa.pdf"),
               beta_plot_list_core_or_not[[core_or_not]][[tax]][[p]][[q]], width = width, height = height)
      }
      
      if (q == "control") { next }
      ggsave(paste0(hist_path,"/beta.", tax,".", p, ".", q, ".hist.pdf"),
             beta_dist[[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      for (core_or_not in names(beta_plot_list_core_or_not)) {
        ggsave(paste0(hist_path_core_or_not,"/beta_core.",core_or_not, ".", tax,".", p, ".", q, ".hist.pdf"),
               beta_dist_core_or_not[[core_or_not]][[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      }
    }
  }
  write_csv(rownames_to_column(beta_dist[[tax]]$dist_df, "Sample_ID"),
            paste0("output/R/beta/beta_all/beta_dist_", tax,".csv"))
  for (core_or_not in names(beta_plot_list_core_or_not)) {
    write_csv(rownames_to_column(beta_dist_core_or_not[[core_or_not]][[tax]]$dist_df, "Sample_ID"),
              paste0("output/R/beta/beta_core_or_not/beta_dist_core_or_not_", core_or_not, ".", tax,".csv"))
  }
}

# Absolute counts
for (tax in names(beta_abs_plot_list)) {
  pcoa_path <- paste0("output/R/beta/beta_all/", tax, "_pcoa")
  pcoa_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_pcoa")
  hist_path <- paste0("output/R/beta/beta_all/", tax, "_hist")
  hist_path_core_or_not <- paste0("output/R/beta/beta_core_or_not/", tax, "_hist")
  system(paste0("mkdir -p ", pcoa_path))
  system(paste0("mkdir -p ", hist_path))
  system(paste0("mkdir -p ", pcoa_path_core_or_not))
  system(paste0("mkdir -p ", hist_path_core_or_not))
  for (p in names(beta_abs_plot_list[[tax]])) {
    if (p=="core") {
      width <- 17
      height <- 15
    } else {
      width <- 17
      height <- 5
    }
    for (q in names(beta_abs_plot_list[[tax]][[p]])) {
      ggsave(paste0(pcoa_path,"/beta_abs.", tax,".", p, ".", q, ".pcoa.pdf"),
             beta_abs_plot_list[[tax]][[p]][[q]], width = width, height = height)
      for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
        ggsave(paste0(pcoa_path_core_or_not,"/beta_abs_core.", core_or_not, ".", tax,".", p, ".", q, ".pcoa.pdf"),
               beta_abs_plot_list_core_or_not[[core_or_not]][[tax]][[p]][[q]], width = width, height = height)
      }
      
      if (q == "control") { next }
      ggsave(paste0(hist_path,"/beta_abs.", tax,".", p, ".", q, ".hist.pdf"),
             beta_abs_dist[[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
        ggsave(paste0(hist_path_core_or_not,"/beta_abs_core.",core_or_not, ".", tax,".", p, ".", q, ".hist.pdf"),
               beta_abs_dist_core_or_not[[core_or_not]][[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
      }
    }
  }
  write_csv(rownames_to_column(beta_abs_dist[[tax]]$dist_df, "Sample_ID"),
            paste0("output/R/beta/beta_all/beta_abs_dist_", tax,".csv"))
  for (core_or_not in names(beta_abs_plot_list_core_or_not)) {
    write_csv(rownames_to_column(beta_abs_dist_core_or_not[[core_or_not]][[tax]]$dist_df, "Sample_ID"),
              paste0("output/R/beta/beta_core_or_not/beta_abs_dist_core_or_not_", core_or_not, ".", tax,".csv"))
  }
}

# ## Heatmaps ####
system("mkdir -p output/R/heatmaps/contigs")
for (tax in names(order_heatmaps_row)) {
  ggsave(paste0("output/R/heatmaps/contigs/", tax, ".pdf"),
         order_heatmaps_row[[tax]],
         width = 20 + plotted_samples[[tax]]/5, height =2 + plotted_contigs[[tax]]/5,
         limitsize = FALSE)
}
system("mkdir -p output/R/heatmaps/tax_collapse")
for (tl in names(tax_collapse_heatmaps)) {
  system(paste0("mkdir -p output/R/heatmaps/tax_collapse/",tl))
  for (tax in names(tax_collapse_heatmaps[[tl]])) {
    tax_clean <- str_replace_all(tax, "\\|", "_")
    ggsave(paste0("output/R/heatmaps/tax_collapse/", tl,"/", tax_clean, ".pdf"),
           tax_collapse_heatmaps[[tl]][[tax]],
           width = 20 + plotted_samples_tax_collapse[[tl]][[tax]]/5, height =2+plotted_contigs_tax_collapse[[tl]][[tax]]/5,
           limitsize = FALSE)
  }
}

## TPM ####
# Hostgroup
system("mkdir -p output/R/relative_abundance/relative_abundance_overall/")
for (p in names(hostgroup_hist)) {
  ggsave(paste0("output/R/relative_abundance/relative_abundance_overall/hostgroup_hist_", p, ".pdf"),
         hostgroup_hist[[p]], width=8, height=4)
}
for (tl in names(average_tpm$Host_group$plots)) {
  ggsave(paste0("output/R/relative_abundance/relative_abundance_overall/average_TPM_Host_groups.",tl,".pdf"),
         average_tpm$Host_group$plots[[tl]], width=15, height=8)
  write_csv(average_tpm$Host_group$tibbles[[tl]],
            paste0("output/R/relative_abundance/relative_abundance_overall/average_TPM_Host_groups.",tl,".csv"))
}

# Core or not
for (core_or_not in names(core_or_not_hist)) {
  ggsave(paste0("output/R/relative_abundance/relative_abundance_overall/averageTPM_core_",core_or_not,".",".pdf"),
         core_or_not_hist[[core_or_not]], width=15, height=8)
}

## Average TPM
system("mkdir -p output/R/relative_abundance/relative_abundance_by_metavar/")
for (tax in names(average_tpm$all)) {
  for (metavar in names(average_tpm$all[[tax]]$plots)) {
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar/relative_abundance_", tax, ".", metavar, ".pdf"),
           average_tpm$all[[tax]]$plots[[metavar]],  width=40, height=15)
    write_csv(average_tpm$all[[tax]]$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_metavar/relative_abundance_", tax, ".", metavar, ".csv"))
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar/HostGroup_relative_abundance_", ".", metavar, ".pdf"),
           average_tpm$Host_group$plots[[metavar]],  width=40, height=15)
    write_csv(average_tpm$Host_group$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_metavar/HostGroup_relative_abundance_", ".", metavar, ".csv"))
  }
}

## Core TPM
system("mkdir -p output/R/relative_abundance/core_TPM/")
for (mvar in names(core_tpm_stats)) {
  ggsave(paste0("output/R/relative_abundance/core_TPM/core_TPM.", mvar, ".pdf"),
         core_tpm_plots[[mvar]], width = 6, height = 6)
  write_csv(core_tpm_stats[[mvar]], paste0("output/R/relative_abundance/core_TPM/core_TPM.", mvar, ".csv"))
}

system("mkdir -p output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/")
for (tax in names(average_tpm_core_or_not$yes)) {
  for (core_or_not in c("yes", "no")) { # <-- This is terrible, but I don't have time to fix it...
    for (metavar in names(average_tpm_core_or_not[[core_or_not]][[tax]]$plots)) {
      wid = 12
      if (metavar == "Sample_ID") {
        wid = 45
      }
      ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/relative_abundance_", core_or_not, ".", tax, ".", metavar, ".pdf"),
             average_tpm_core_or_not[[core_or_not]][[tax]]$plots[[metavar]],  width=wid, height=5)
      write_csv(average_tpm_core_or_not[[core_or_not]][[tax]]$tibbles[[metavar]],
                paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/relative_abundance_", core_or_not, ".", tax, ".", metavar, ".csv"))
      ggsave(paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/Core_or_not_relative_abundance.", metavar, ".pdf"),
             average_tpm_core_or_not$Core_or_not$plots[[metavar]],  width=wid, height=5)
      write_csv(average_tpm_core_or_not$Core_or_not$tibbles[[metavar]],
                paste0("output/R/relative_abundance/relative_abundance_by_metavar_core_or_not/Core_or_not_relative_abundance.", metavar, ".csv"))
    }
  }
}

## Prevalence ####
system("mkdir -p output/R/prevalence/")
for (tl in names(prevalence_histo)) {
  if (tl=="Bee_pools") {
    wid=20
    plot <- prevalence_histo[[tl]]$plot # + scale_x_continuous(breaks = seq(3,150,3))
  } else if (tl=="Hives") {
    wid=12
    plot <- prevalence_histo[[tl]]$plot 
  } else {
    wid=8
    plot <- prevalence_histo[[tl]]$plot
  }
  ggsave(paste0("output/R/prevalence/prevalence.",tl,".pdf"),
         plot, width=wid, height=5)
  write_csv(prevalence_histo[[tl]]$table, paste0("output/R/prevalence/prevalence.",tl,".csv"))
}

## TPM per country ####
system("mkdir -p output/R/relative_abundance/relative_abundance_by_country")
for (core_or_not in names(average_tpm_per_country)) {
  for (tl in names(average_tpm_per_country[[core_or_not]])) {
    for (country in names(average_tpm_per_country[[core_or_not]][[tl]])) {
      ggsave(paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_core.", core_or_not, ".", tl, ".", country, ".pdf"),
             average_tpm_per_country[[core_or_not]][[tl]][[country]], height = 15, width = 35)
    }
  }
}
for (country in names(average_tpm_per_country_metavar)) {
  for (metavar in names(average_tpm_per_country_metavar[[country]]$plots)) {
    ggsave(paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_metavar.", metavar, ".", country, ".pdf"),
           average_tpm_per_country_metavar[[country]]$plots[[metavar]], height = 5, width = 15)
    write_csv(average_tpm_per_country_metavar[[country]]$tibbles[[metavar]],
              paste0("output/R/relative_abundance/relative_abundance_by_country/relative_abundance_metavar.", metavar, ".", country, ".csv"))
  }
}

## Prevalence per country ####
system("mkdir -p output/R/prevalence/prevalance_by_country")
for (core_or_not in names(prevalence_plots_per_country)) {
  for (tl in names(prevalence_plots_per_country[[core_or_not]])) {
    for (country in names(prevalence_plots_per_country[[core_or_not]][[tl]])) {
      ggsave(paste0("output/R/prevalence/prevalance_by_country/prevalence_core.", core_or_not, ".", tl, ".", country, ".pdf"),
             prevalence_plots_per_country[[core_or_not]][[tl]][[country]],
             width = 30, height = 20)
    }
  }
}

## Upset ####
pdf(file="output/R/venns/country_upset.pdf",
    width = 25, height = 5) 
country_upset_plot
dev.off()


## Genera per family
system("mkdir -p output/R/gener_per_family/")
ggsave("output/R/gener_per_family/genera_per_family.pdf", genera_per_family_plot,
       width = 6, height = 6)
write_csv(genera_per_family_stats, "output/R/gener_per_family/genera_per_family.csv")


end_time <- Sys.time()
end_time - start_time

