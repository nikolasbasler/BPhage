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

bphage_microvirus_contigs <- read_lines("data/bphage.microviridae.contigs") # Where is this generated?

bphage_microvirus_taxonomy <- read.csv("data/bphage.microvirus.taxonomy.csv", sep=";") # Where is this generated?

classification_gnmd <- read.csv("output/R/phage.filt.gnmd.classification.csv") %>%
  mutate(contig_length = contig_length/1000) %>%
  rename(length_kb = contig_length)

present_in_all_countries <- read_lines("data/core_contigs.txt") # Where is this generated?

prevalences_tables <- list()
for (thing in c("Bee_pools", "Hives", "Countries")) {
  prevalences_tables[[thing]] <- read.csv(paste0("data/prevalence_tables/prevalence.",thing,".csv")) %>%
    select(contig, prevalence_abs) %>%
    rename(!!paste0("Prevalence_",thing) := prevalence_abs) %>% 
    tibble()
}

host_group_df <- read.csv("data/host_groups.csv") # This df is generated in host.R and saved into data/ for convenience, so it can be used here already.

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
  left_join(., bphage_microvirus_taxonomy, by ="contig") %>% 
  mutate(Subgroup = ifelse(Class != "Faserviricetes|Huolimaviricetes|Malgrandaviricetes", NA, Subgroup)) %>% # Only call microvirus what was classified in the appropriate class by vcontact
  mutate(Family = ifelse(!is.na(Subgroup), paste0("Microvirus_", Subgroup), Family)) %>%
  mutate(Order = ifelse(contig %in% bphage_microvirus_contigs & Class == "Faserviricetes|Huolimaviricetes|Malgrandaviricetes", "Microviruses", Order)) %>% # Only call microvirus what was classified in the appropriate class by vcontact
  select(-Subgroup) %>%
  mutate(lowest_taxon = apply(., 1, function(row) tail(row[row != ""], 1))) %>%
  mutate(across(everything(), ~ifelse(. == "", "Unclassified", .))) %>%
  inner_join(., classification_gnmd[c("contig", "provirus", "proviral_length", 
                                      "gene_count", "viral_genes", "host_genes",
                                      "checkv_quality", "miuvig_quality", 
                                      "completeness", "completeness_method",
                                      "contamination", "kmer_freq", "warnings")],
             by = "contig") %>%
  mutate(predicted_genome_length = length_kb / (completeness/100)) %>%
  left_join(., host_group_df, by = "contig") %>%
  mutate(Core = ifelse(contig %in% present_in_all_countries, "yes", "no")) %>%
  mutate(Order_group = Order) %>%
  mutate(Order_group = ifelse(str_detect(Order_group, "novel_order.*of_Caudoviricetes"), "Novel_Caudoviricetes_order", Order_group)) %>%
  mutate(Order_group = ifelse(str_detect(Order_group, "novel_order.*of_Tokiviricetes"), "Novel_Tokiviricetes_order", Order_group)) %>%
  mutate(Family_group = Family) %>%
  mutate(Family_group = ifelse(str_detect(Family_group, "novel_family.*of_Caudoviricetes"), "Novel_Caudoviricetes_family", Family_group)) %>%
  mutate(Family_group = ifelse(str_detect(Family_group, "novel_family.*of_Tokiviricetes"), "Novel_Tokiviricetes_family", Family_group)) %>%    
  mutate(Family_group = ifelse(!str_detect(Family_group, "Novel") & !str_detect(Family_group, "Micro") & Family_group !="Unclassified", "ICTV-named", Family_group)) %>%
  mutate(Family_group = ifelse(str_detect(Family_group, "Micro"), "Microvirus_family", Family_group)) %>%
  mutate(Family_group = ifelse(Family_group == "Unclassified" & Order == "Microviruses", "Unclassified_Microvirus", Family_group)) %>%
  mutate(Family_group = ifelse(Family_group == "Unclassified" & Order != "Microviruses", "Other_unclassified", Family_group)) %>%

  left_join(., prevalences_tables$Bee_pools, by = "contig") %>%
  left_join(., prevalences_tables$Hives, by = "contig") %>%
  left_join(., prevalences_tables$Countries, by = "contig") 

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
# Box plots of absolute viral loads ####
vlp_stats <- metadata %>%
  select(VLPs_per_ul) %>%
  filter(!is.na(VLPs_per_ul)) %>%
  summarise(q25 = round(quantile(VLPs_per_ul)[2]),
            median = round(median(VLPs_per_ul)),
            q75 = round(quantile(VLPs_per_ul)[4]),
            sd = round(sd(VLPs_per_ul)),
            norm_sd = round(sd(VLPs_per_ul) / mean(VLPs_per_ul), digits = 3)) %>%
  mutate(q25 = paste0("q25: ", as.character(q25)),
         median = paste0("median: ", as.character(median)),
         q75 = paste0("q75: ", as.character(q75)),
         sd = paste0("sd: ", as.character(sd)),
         norm_sd = paste0("norm_sd: ", as.character(norm_sd))) %>%
  unlist()

vlp_box <- metadata %>%
  select(VLPs_per_ul) %>%
  filter(!is.na(VLPs_per_ul)) %>%
  ggplot(aes(y=VLPs_per_ul)) +
  geom_boxplot() +
  annotate("text", x=0.4, y=seq(400000, 800000, 100000), label = vlp_stats, adj = "right") 

vlp_hist <- metadata %>%
  select(VLPs_per_ul) %>%
  filter(!is.na(VLPs_per_ul)) %>%
  ggplot(aes(x=VLPs_per_ul)) + 
  geom_histogram() +
  annotate("text", x=600000, y=seq(10, 30, 5), label = vlp_stats, adj = "right")

vlp_overview <- wrap_plots(list(vlp_box, vlp_hist))

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
taxlevels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", 
               "Realm", "Host_group", "Core", "Order_group", "Family_group",
               "Prevalence_Bee_pools", "Prevalence_Hives", "Prevalence_Countries")
phage_ab <- list()
phage_ab <- taxlevels %>%
  set_names() %>%
  map(~tax_sum(., ab_table=phage_abundance,
               classif = classification))

phage_ab_hostg_core <- list()
phage_ab_hostg_core$all <- phage_ab$Host_group # This is redundant but convenient
for (core_or_not in c("yes", "no")) {
  classif_filt <- classification %>%
    filter(Core == core_or_not)
  phage_ab_hostg_core[[core_or_not]] <- tax_sum("Host_group", ab_table=phage_abundance,
                 classif = classif_filt)
}

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

phage_tpm_hostg_core <- list()
for (group in names(phage_ab_hostg_core)) {
  phage_tpm_hostg_core[[group]] <- calc_tpm(abtable = phage_ab_hostg_core[[group]], 
                                                  level = "Host_group", 
                                                  lengths_df = phage_lengths$Host_group)
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
# Pretty pies # Out-sourced because of a lot of hard-coding

source("scripts/R/helpers/pretty_pies.R") 

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

completeness <- list()
# median_completeness <- unique(median(classification$completeness, na.rm = TRUE))
completeness$plot  <- classification %>%
  filter(completeness >= 50) %>%
  ggplot(aes(x = completeness)) +
  geom_histogram(binwidth = 2) +
  # geom_vline(xintercept = median_completeness, linetype = "dashed") +
  labs(x = "Completeness (%)", y = "Genome count") +
  theme_pubr()

completeness$tibble <- classification %>%
  mutate(completeness_category = case_when(
    completeness >= 100 ~ "complete",
    completeness >= 90  ~ "near_complete_(90%)",
    completeness >= 50  ~ "partial",
    TRUE ~ "NA")) %>%
  select(completeness, completeness_category) %>%
  group_by(completeness_category) %>%
  summarise(count = n(), 
            proportion = n() / 2343)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Set min seq thresholds for rarefaction. ####

# 50% horizontal coverage cutoff
# count_stats <- report_stats(df = phage_abundance,
#                             thresholds=c(2689, 1002)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
# min_seq_count <- 2713 # For quantile_90_ratio < 1000. 11 discarded samples (10 mid, 1 ile).

# 70% horizontal coverage cutoff
count_stats <- report_stats(df = phage_abundance,
                            thresholds=c(2437, 1548, 775)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
count_stats
# min_seq_count <- 2442 # For quantile_90_ratio < 1000. 12 discarded samples (9 mid, 3 ile).
min_seq_count <- 1548 # Looks like a natural breaking point. 4 discarded samples (all mid).
discarded <- discards(count_stats$ratios, min_seq_count)$discarded
lost_bees <- discards(count_stats$ratios, min_seq_count)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

min_seq_count_core_or_not <- list()
count_stats_core_or_not <- list()
count_stats_core_or_not$no <- report_stats(df = phage_ab_core_or_not$no$contig,
                            thresholds=c(961, 1694, 7356))
count_stats_core_or_not$no
min_seq_count_core_or_not$no <- 1694 # Looks like a natural breaking point. 7 discarded samples (5 mid, 2 ile).
discarded <- discards(count_stats_core_or_not$no$ratios, min_seq_count_core_or_not$no)$discarded
lost_bees <- discards(count_stats_core_or_not$no$ratios, min_seq_count_core_or_not$no)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

count_stats_core_or_not$yes <- report_stats(df = phage_ab_core_or_not$yes$contig,
                                     thresholds=c(462, 1073, 2375))
count_stats_core_or_not$yes
min_seq_count_core_or_not$yes <- 1073 # Looks like a natural breaking point. 1 discarded ile sample.
discarded <- discards(count_stats_core_or_not$yes$ratios, min_seq_count_core_or_not$yes)$discarded
lost_bees <- discards(count_stats_core_or_not$yes$ratios, min_seq_count_core_or_not$yes)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Alpha and beta diversity # Out-sourced to easily deactivate for testing, 
# because this part takes by far the longest.

iterations <- 1000
# source("scripts/R/helpers/alpha_beta_rarefaction.R") 

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
# Average TPM ####
met_v <- c("Country", "Season", "Gut_part", "Health", "Sample_ID")
tax_levels <- c("contig", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
average_tpm <- list()
for (tlvl in tax_levels) {
  average_tpm[[tlvl]] <- average_tpm_bar_plot(
    tpm_table = phage_tpm[[tlvl]],
    tl = tlvl,
    title_prefix = "Host group",
    hg = "all",
    meta_vars = met_v)
}

average_tpm_core_or_not <- list()
average_tpm_core_or_not$Core_or_not <- average_tpm_bar_plot(tpm_table = phage_tpm$Core,
                                                            tl = "Core_or_not",
                                                            hg = "Core_or_not",
                                                            meta_vars = met_v,
                                                            hg_or_core = "Core?",
                                                            threshold_for_other = 0)

for (meta_merge in c("Prevalence_Bee_pools", "Prevalence_Hives", "Prevalence_Countries")) {
  average_tpm_core_or_not[[meta_merge]] <- average_tpm_bar_plot(tpm_table = phage_tpm[[meta_merge]],
                                                              tl = "Prevalence",
                                                              hg = meta_merge,
                                                              meta_vars = met_v,
                                                              hg_or_core = "",
                                                              threshold_for_other = 0)
}

average_tpm_core_or_not_taxes <- list()
for (core_or_not in names(phage_tpm_core_or_not)) {
  for (tlvl in tax_levels) {
    average_tpm_core_or_not_taxes[[core_or_not]][[tlvl]] <- average_tpm_bar_plot(
      tpm_table = phage_tpm_core_or_not[[core_or_not]][[tlvl]],
      tl = tlvl,
      hg = core_or_not,
      meta_vars = met_v,
      threshold_for_other=0.01,
      hg_or_core = "Core?")
  }
}
# average_tpm_host_group <- list()
# average_tpm_host_group <- average_tpm_bar_plot(tpm_table = phage_tpm$Host_group,
#                                                tl = "Host_group", 
#                                                hg = "Host_group",
#                                                meta_vars = met_v,
#                                                threshold_for_other = 0,
#                                                hg_or_core = "Host_group")

average_tpm_host_group <- list()
for (hostg in names(phage_tpm_hostg_core)) {
  average_tpm_host_group[[hostg]] <- average_tpm_bar_plot(tpm_table = phage_tpm_hostg_core[[hostg]],
                                                 tl = "Host_group", 
                                                 hg = paste0("Core: ", hostg),
                                                 meta_vars = met_v,
                                                 threshold_for_other = 0,
                                                 hg_or_core = "Host_group")
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
    left_join(., classification, by = "contig") %>%
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
# Prevalenve Venn diagrams # Out-sourced for easy deactivation, because output 
# files are automatically written by the tool.

# source("scripts/R/helpers/venn_run.R") 

here_time <- Sys.time()
here_time - start_time

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Save files # Out-sourced for easy deactivation.

# source("scripts/R/helpers/save_files.R") 


end_time <- Sys.time()
end_time - start_time

