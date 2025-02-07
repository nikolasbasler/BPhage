library(tidyverse)
library(ggpubr)
library(multcompView)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

# phold_all_cds_functions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_all_cds_functions_long_names.tsv.gz")
phold_per_cds_predictions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_per_cds_predictions_long_names.tsv.gz")
extended_phold_per_cds_predictions <- read.delim("output/core_contig_refinement/extended_contigs_phold/phold_per_cds_predictions_long_names.tsv") %>%
  mutate(
    contig_id = str_extract(contig_id, paste(rep("([^_]+)", 11), collapse = "_")),
    cds = sapply(str_split(cds_id, "_"), function(x) paste(tail(x, 2), collapse = "_"))
  ) %>%
  mutate(cds_id = paste(contig_id, cds, sep = "_")) %>%
  select(-cds)

phold_predictions_with_extensions <-  phold_per_cds_predictions %>% 
  filter(!contig_id %in% unique(extended_phold_per_cds_predictions$contig_id)) %>%
  bind_rows(., extended_phold_per_cds_predictions)



### Boxplots with differences between core and non-core genes ####

# SOMEHOW MAKE IT SO THAT ONLY BOXES WITH SIGNIFICANT PWCs ARE DRAWN 
# I hard-coded this now. Automate it.

# test_metrics <- c("function.", "product", "transl_table")
# annotation_core_plots <- list()
# for (metric in test_metrics) {
#   annotation_core_plots[[metric]] <- phold_predictions_with_extensions %>%
#     mutate(transl_table = as.character(transl_table)) %>%
#     filter(str_starts(contig_id, "NODE")) %>%
#     group_by(contig_id, .data[[metric]]) %>%
#     count() %>%
#     ungroup() %>% 
#     left_join(., classification, by = join_by(contig_id == contig)) %>%
#     filter(!is.na(Core)) %>%
#     group_by(Core) %>%
#     mutate(genes_per_contig_kb = n / length_kb) %>%
#     mutate(genes_per_genome_kb = n / predicted_genome_length) %>%
#     ungroup() %>%
#     select(all_of(metric), genes_per_contig_kb, genes_per_genome_kb, Core) %>%
#     pivot_longer(-c(all_of(metric), Core)) %>% 
#     ggplot(aes(x = all_of(metric), y = value, fill = Core)) +
#     geom_boxplot(aes(x = .data[[metric]], y = value)) +
#     facet_wrap(~name) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = metric)
# }
# annotation_core_plots

# Pairwise comparisons between products doesn't work, probably because there is 
# variance of 0 between the groups.
# test_metrics <- c("function.", "product", "transl_table")
# annotation_core_plots <- list()
# for (metric in test_metrics) {
#   plot_tbl <- phold_predictions_with_extensions %>%
#     mutate(transl_table = as.character(transl_table),
#            function. = ifelse(function. == "", "unknown function", function.)) %>%
#     filter(str_starts(contig_id, "NODE")) %>%
#     group_by(contig_id, .data[[metric]]) %>%
#     count() %>%
#     ungroup() %>%
#     left_join(., classification, by = join_by(contig_id == contig)) %>%
#     filter(Class == "Caudoviricetes",
#            completeness >= 90) %>%
#     group_by(.data[[metric]], Core) %>%
#     filter(sum(n) > 1) %>%
#     group_by(.data[[metric]]) %>%
#     filter(any(Core == "yes") & any(Core == "no")) %>%
#     ungroup() # %>%
#     # filter(.data[[metric]] != "hypothetical protein")
#   
#   # plot_tbl %>%
#   #   group_by(.data[[metric]], Core) %>%
#   #   summarize(variance = var(n), count = n(), .groups = "drop") %>%
#   #   filter(count < 2 | variance == 0)
# 
#    annotation_core_plots[[metric]] <- plot_tbl %>%
#      ggplot(aes(x = .data[[metric]], y = n, fill = Core)) +
#      geom_boxplot(aes(x = .data[[metric]], y = n)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = metric) +
#      geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE, group.by = "x.var")
   
   # if (metric == "product") {
   #   significant_products <- c("anti-repressor Ant", "baseplate hub", "baseplate spike",
   #                             "baseplate wedge subunit", "beta-propeller repeat protein",
   #                             "DefenseFinder protein", "DNA binding protein", "DNA polymerase",
   #                             "GTP-binding domain", "head decoration",
   #                             "head maturation protease", "head morphogenesis",
   #                             "head scaffolding protein", "head-tail adaptor", "holin",
   #                             "homing endonuclease", "integrase", "lipoprotein", "major head protein",
   #                             "Mu Gam-like end protection", "nucleotide kinase",
   #                             "portal protein", "RecT-like ssDNA annealing protein",
   #                             "replication initiation protein", "RIIB lysis inhibitor",
   #                             "Rz-like spanin", "single strand DNA binding protein",
   #                             "structural protein", "tail assembly chaperone",
   #                             "tail chaperone protein", "tail completion or Neck1 protein",
   #                             "tail fiber protein", "tail length tape measure protein",
   #                             "tail protein", "tail sheath", "tail terminator",
   #                             "terminase large subunit", "terminase small subunit",
   #                             "transcriptional repressor", "virion structural protein")
   #   plot_tbl_filt <- plot_tbl %>%
   #     filter(product %in% significant_products)
   #   plot_tbl_filt %>% select(product) %>% distinct() %>% unlist(use.names = FALSE) %>% sort()
   # 
   #   annotation_core_plots$significant_products <- plot_tbl_filt %>%
   #     ggplot(aes(x = .data[[metric]], y = n, fill = Core)) +
   #     # ggplot(aes(x = all_of(metric), y = n)) +
   #     geom_boxplot(aes(x = .data[[metric]], y = n)) +
   #     # facet_wrap(~name) +
   #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   #     labs(title = metric) +
   #     geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE, group.by = "x.var")
# 
#    }
# }


### Bar chart for core PHROGs ####
moron_bar_colors <- rev(c("lightgrey", "#ef8f01", "#555555"))

contigs_per_phrog <- phold_predictions_with_extensions %>%
  filter(contig_id %in% present_in_all_countries) %>%
  group_by(function.) %>%
  summarise(number_of_contigs = n_distinct(contig_id)) %>%
  arrange(desc(number_of_contigs))

phrog_tibble <- phold_predictions_with_extensions %>%
  filter(contig_id %in% present_in_all_countries) %>%
  group_by(function.) %>%
  summarise(count = n()) %>%
  mutate(collapsed_cat = "Other PHROG\ncategories") %>%
  mutate(collapsed_cat = ifelse(function. == "unknown function", "unknown function", collapsed_cat),
         collapsed_cat = ifelse(function. == "moron, auxiliary metabolic gene and host takeover", '"moron", AMG\nand host takeover', collapsed_cat)
  ) %>%
  group_by(collapsed_cat) %>%
  summarise(gene_count = sum(count)) %>%
  mutate(collapsed_cat = factor(collapsed_cat, levels = c("Other PHROG\ncategories", '"moron", AMG\nand host takeover', "unknown function")))

phrog_bar <- phrog_tibble %>%
  ggplot(aes(x = "", y = gene_count, fill = collapsed_cat)) +
  geom_bar(stat = "identity", color = "black") +
  theme_void() +
  theme(
    legend.text = element_text(hjust = 1),
    legend.margin=margin(0,20,10,0)) +
  scale_fill_manual(values = moron_bar_colors) +
  labs(fill = "PHROG category") +
  theme(legend.spacing.x = unit(0.25, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))

phrog_bar_horizontal <- phrog_bar + 
  coord_flip() +
  theme(legend.position = "top",
        legend.text = element_text(hjust = 0)) +
  guides(fill = guide_legend(reverse=T))

moron_with_groups <- phold_predictions_with_extensions %>%
  filter(contig_id %in% present_in_all_countries) %>%
  filter(function. == "moron, auxiliary metabolic gene and host takeover") %>% 
  # group_by(product) %>%
  # summarise(count = n()) %>%
  # arrange(desc(count)) %>%
  mutate(group = NA) %>%
  mutate(group = ifelse(product == "membrane protein" |
                          product == "PAAR motif of membran proteins" |
                          product == "ABC transporter" |
                          product == "membrane associated protein",
                        "Membrane protein", group),
         group = ifelse(product == "phosphoadenosine phosphosulfate reductase", "Sulfur metabolism", group),
         group = ifelse(product == "gam-like host nuclease inhibitor" |
                          product == "anti-restriction protein" |
                          product == "anti-CRISPR protein" |
                          product == "host nuclease inhibitor" |
                          product == "DarB-like antirestriction" |
                          product == "Lar-like restriction alleviation protein" |
                          product == "anti-sigma factor",
                        "Host takeover", group),
         group = ifelse(product == "Doc-like toxin"|
                          product == "MazF-like growth inhibitor" |
                          product == "toxin" |
                          product == "toxin-antitoxin system HicB-like" |
                          product == "RelE-like toxin" |
                          product == "plasmid antitoxin with HTH domain" |
                          product == "ribonuclease toxin of AT system",
                        "Toxin-antitoxin system", group),
         group = ifelse(product == "abortive infection resistance protein" |
                          product == "superinfection exclusion" |
                          product == "SieB superinfection exclusion" |
                          product == "DefenseFinder protein" |
                          product == "superinfection exclusion Sie-like",
                        "Superinfection exclusion", group),
         group = ifelse(product == "beta-lactamase-inhibitor protein BLIP" |
                          product == "tellurite resistance",
                        "Antibiotic resistance", group),
         group = ifelse(product == "queuine tRNA-ribosyltransferase" |
                          product == "levanase" |
                          product == "PnuC-like nicotinamide mononucleotide transport" |
                          product == "VFDB virulence factor protein",
                        "other", group)
  )

contigs_with_TA <- moron_with_groups %>%
  filter(group == "Toxin-antitoxin system") %>%
  group_by(contig_id) %>%
  summarise(TA_per_contig = n()) %>%
  left_join(., classification, by = join_by(contig_id == contig)) %>%
  select(contig_id, TA_per_contig, length_kb, lowest_taxon, completeness, Host_group) %>%
  arrange(desc(TA_per_contig))

moron_grouping <- moron_with_groups %>%
  select(product, group) %>%
  group_by(product) %>%
  mutate(count = n()) %>%
  arrange(desc(count)) %>%
  distinct()

moron_tibble <- moron_grouping %>%
  group_by(group) %>%
  summarise(group_count = sum(count)) %>%
  arrange(group_count) %>%
  mutate(group = factor(group, levels = c("other", "Antibiotic resistance", "Sulfur metabolism", 
                                              "Toxin-antitoxin system", "Host takeover", 
                                              "Membrane protein","Superinfection exclusion")))

genes_in_core <- phold_predictions_with_extensions %>%
  filter(contig_id %in% present_in_all_countries) %>%
  group_by(product) %>%
  summarise(gene_count = n()) %>%
  arrange(desc(gene_count))

### Pie chart for core morons ####
moron_pie_colors <- c("#555555", "#FFDAB9", "#FFA07A", "#FFC300", "#D2691E", "#8B4513", "#5E280C")
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


## Sulfur metabolism gene

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv")

sulfur_phage_annot <- phold_predictions_with_extensions %>%
  # filter(str_detect(product, "sulf")) %>%
  filter(product == "phosphoadenosine phosphosulfate reductase") %>% 
  inner_join(classification, ., by = join_by("contig" == "contig_id"))

sulf_meta <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "TPM") %>%
  mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE)) %>%
  left_join(., metadata, by = "Sample_ID") %>% 
  select(contig, Sample_ID, TPM, sulf_metabolism, Hive_ID, Country, Season, Gut_part, Health) %>%
  filter(Gut_part == "rec") %>% # Focus on rectum for more stable values (some midgut samples have very few phages)
  group_by(Sample_ID, sulf_metabolism) %>% 
  mutate(sulf_tpm = sum(TPM)) %>%
  ungroup()

sulf_positive_hives <- sulf_meta %>%
  filter(sulf_tpm > 0) %>%
  group_by(Country) %>%
  summarise(hive_count = n_distinct(Hive_ID),
            sulf_positive_hives = n_distinct(Hive_ID[sulf_metabolism])) %>%
  mutate(sulf_positive_prop = sulf_positive_hives / hive_count)

meta_variables <- c("Hive_ID", "Country", "Season", "Health")
sulf_stats <-  list()
tpm_KW <- list()
tpm_pwc <- list()
sulf_tpm_plots <- list()
genomes_KW <- list()
genomes_pwc <- list()
sulf_genome_prop_plots <- list()
for (met_v in meta_variables) {
  tpm <- sulf_meta %>%
    filter(sulf_metabolism) %>%
    group_by(.data[[met_v]]) %>%
    mutate(mean_sulf_tpm = mean(sulf_tpm)) %>%
    filter(mean_sulf_tpm > 0) %>%
    ungroup() %>%
    select(Sample_ID, all_of(meta_variables), sulf_tpm, mean_sulf_tpm) %>%
    distinct() %>%
    arrange(desc(mean_sulf_tpm))
  tpm_KW[[met_v]] <- kruskal.test(tpm$sulf_tpm ~ tpm[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))
  # See if I can make this work. Test with Season to see the problem. Letters aren't assigned to every group.
  # pwc_p <- pairwise.wilcox.test(tpm$sulf_tpm, tpm[[met_v]], p.adjust.method = "BH")
  # all_groups <- unique(c(rownames(pwc_p$p.value), colnames(pwc_p$p.value)))
  # full_p_matrix <- as.data.frame(matrix(1, nrow = length(all_groups), ncol = length(all_groups), 
  #                                       dimnames = list(all_groups, all_groups)))
  # full_p_matrix[rownames(pwc_p$p.value), colnames(pwc_p$p.value)] <- pwc_p$p.value
  # multcompLetters(full_p_matrix)$Letters
  # tpm_pwc[[met_v]] <- multcompLetters(pwc_p$p.value)$Letters

  genomes <- sulf_meta %>%
    filter(TPM >0) %>%
    group_by(Sample_ID, sulf_metabolism) %>% 
    mutate(sulf_genomes = n()) %>%
    group_by(Sample_ID) %>%
    mutate(sample_sulf_genome_props = sulf_genomes / n ()) %>% 
    group_by(.data[[met_v]]) %>%
    mutate(sulf_genome_count = sum(sulf_metabolism),
           sulf_genome_prop = mean(sulf_metabolism)) %>%
    ungroup() %>%
    filter(sulf_metabolism) %>%
    select(Sample_ID, all_of(meta_variables), sample_sulf_genome_props, sulf_genome_prop) %>%
    distinct() %>%
    arrange(desc(sulf_genome_prop))
  genomes_KW[[met_v]] <- kruskal.test(genomes$sample_sulf_genome_props ~ genomes[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))
  pwc_p <- pairwise.wilcox.test(genomes$sample_sulf_genome_props, genomes[[met_v]], p.adjust.method = "BH")
  genomes_pwc[[met_v]] <- multcompLetters(pwc_p$p.value)$Letters
  
  sulf_tpm_plots[[met_v]] <- tpm %>%
    # mutate(!!met_v := factor(.data[[met_v]], levels = unique(tpm[[met_v]]))) %>%
    ggplot(aes(x = .data[[met_v]], y = sulf_tpm)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25) +
    # scale_y_continuous(trans='log10') + # This would change the picture because 0s would be ignored.
    ggtitle(paste0("KW p-value: ", tpm_KW[[met_v]]$p.value)) # +
    # scale_y_continuous(limits = c(0,0.1))

  sulf_genome_prop_plots[[met_v]] <- genomes %>%
    # mutate(!!met_v := factor(.data[[met_v]], levels = unique(genomes[[met_v]]))) %>%
    ggplot(aes(x = .data[[met_v]], y = sample_sulf_genome_props)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25) +
    ggtitle(paste0("KW p-value: ", genomes_KW[[met_v]]$p.value))
  
  if (met_v %in% c("Hive_ID", "Season")) {
    sulf_tpm_plots[[paste0(met_v,"_facet")]] <- sulf_tpm_plots[[met_v]] +
      facet_wrap(~Country, scales = "free_x") +
      ggtitle("")
    sulf_genome_prop_plots[[paste0(met_v,"_facet")]] <- sulf_genome_prop_plots[[met_v]]  +
      facet_wrap(~Country, scales = "free_x") +
      ggtitle("")
  }
  
  summarised_tpm <- tpm %>%
    select(all_of(met_v), mean_sulf_tpm) %>%
    distinct()
  summarised_genomes <- genomes %>% 
    select(all_of(met_v), sulf_genome_prop) %>% 
    distinct()
  sulf_stats[[met_v]] <- left_join(summarised_tpm, summarised_genomes, by = met_v)
}



### Save files ####
system("mkdir -p output/R/gene_content/sulfur")
# for (plot in names(annotation_core_plots)) {
#   wid <- 12
#   hig <- 8
#   if (plot == "product") {
#     wid <- 50
#     hig <- 10
#   }
#   ggsave(paste0("output/R/gene_content/gene_content.", plot, ".pdf"), 
#          annotation_core_plots[[plot]], width = wid, height = hig, limitsize = FALSE)
# }

ggsave("output/R/gene_content/core_phrog_bar.pdf",
       phrog_bar, width = 2.5, height = 6)
ggsave("output/R/gene_content/core_phrog_bar_horiziontal.pdf",
       phrog_bar_horizontal, width = 5.5, height = 1)
write_csv(phrog_tibble, "output/R/gene_content/core_phrog_bar.csv")
write_csv(contigs_per_phrog, "output/R/gene_content/core_contigs_per_phrog.csv")

ggsave("output/R/gene_content/core_moron_pie.pdf",
       moron_pie,  width = 4.5, height = 4.5)
write_csv(moron_grouping, "output/R/gene_content/core_moron_pie.csv")
write_csv(contigs_with_TA, "output/R/gene_content/core_TA_contigs.csv")
write_csv(genes_in_core, "output/R/gene_content/core_all_products.csv")

## Sulfur

write_delim(sulf_positive_hives, "output/R/gene_content/sulfur/sulfur_positive_hives.tsv", delim = "\t ")
for (meta in names(sulf_stats)) {
  write_delim(sulf_stats[[meta]], paste0("output/R/gene_content/sulfur/sulfur_stats.", meta, ".tsv"), delim = "\t")
  write_delim(tpm_KW[[meta]], paste0("output/R/gene_content/sulfur/tpm_krusk.", meta, ".tsv"), delim = "\t")
  write_delim(genomes_KW[[meta]], paste0("output/R/gene_content/sulfur/genome_count_krusk.", meta, ".tsv"), delim = "\t")
}
for (meta in names(sulf_tpm_plots)) {
  ggsave(paste0("output/R/gene_content/sulfur/tpm.", sulf_tpm_plots[[meta]], ".pdf"),
         width = 7, height = 7)
  ggsave(paste0("output/R/gene_content/sulfur/genome_count.", sulf_genome_prop_plots[[meta]], ".pdf"),
         width = 7, height = 7)
  
  as.matrix() %>% t() %>% as_tibble()
}

