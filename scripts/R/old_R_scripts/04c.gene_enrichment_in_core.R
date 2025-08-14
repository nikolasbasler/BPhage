library(tidyverse)
library(ggrepel)

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

core_genome_count <- length(present_in_all_countries)
# noncore_genome_count <- nrow(classification) - core_genome_count
caudo_genomes <- classification %>% 
  filter(Class == "Caudoviricetes") %>%
  select(contig) %>%
  unlist(use.names = FALSE) %>%
  as.vector()

noncore_caudo_count <- setdiff(caudo_genomes, present_in_all_countries) %>% 
  length()

gene_counts_by_core <- phold_predictions_with_extensions %>%
  filter(contig_id %in% caudo_genomes) %>% 
  mutate(Core = ifelse(contig_id %in% present_in_all_countries, "yes", "no")) %>%
  select(contig_id, product, Core) %>% 
  distinct() %>%
  count(product, Core, name = "n_present") %>%
  complete(product, Core = c("yes", "no"), fill = list(n_present = 0)) %>%
  mutate(n_absent = ifelse(Core == "yes", core_genome_count - n_present,
                           noncore_caudo_count - n_present))

contingency_matrices <- list()
ftests <- list()
for (gene in unique(gene_counts_by_core$product)) {
  cont_matr <- gene_counts_by_core %>%
    filter(product == gene) %>%
    select(-product) %>% 
    arrange(desc(Core)) %>%
    column_to_rownames("Core") %>%
    as.matrix()

  contingency_matrices[[gene]] <- cont_matr
  ftests[[gene]] <- fisher.test(cont_matr)
}

# In this arrangement:
# Odds ratio >1: Gene enriched in core
# Odds ratio <1: Gene depleted in core

ftest_summary <- data.frame(p_value = map_dbl(ftests, "p.value"), 
                            odds_ratio = map_dbl(ftests, "estimate"),
                            CI_lower = map_dbl(ftests, ~ .x$conf.int[1]),
                            CI_upper = map_dbl(ftests, ~ .x$conf.int[2])) %>%
  rownames_to_column("gene") %>%
  tibble() %>%
  mutate(p_value = p.adjust(p_value, "BH")) %>%
  rename(adjusted_p = p_value)
# 
# ftest_summary %>%
#   filter(adjusted_p <= 0.05) %>% View()

list_of_significants <- ftest_summary %>%
  filter(adjusted_p <= 0.05) %>%
  ggplot(aes(x = reorder(gene, odds_ratio), y = odds_ratio)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    geom_hline(yintercept = 1, color = "red") +
    coord_flip() +
  scale_y_continuous(trans='log10') +
    labs(x = "Gene", 
         y = "Odds Ratio (with 95% CI)",
         # y = "Odds Ratio",
         title = "Gene Enrichment in Core Viruses") +
    theme_minimal()

# Create the volcano plot
volcano <- ftest_summary %>%
  mutate(log2_OR = log2(odds_ratio),
         neg_log10_p = -log10(adjusted_p)) %>%
  ggplot(aes(x = log2_OR, y = neg_log10_p, color = adjusted_p <= 0.05)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    labs(x = "Log2 Odds Ratio", 
         y = "-Log10(P-value)",
         title = "Gene Enrichment in Core Viruses") +
    theme_minimal() +
    scale_color_manual(values = c("black", "red"), guide = "none")

# Save files
write_delim(ftest_summary, "output/R/gene_content/enrichment_ftests.tsv", delim = "\t")
ggsave("output/R/gene_content/enrichment_significants.pdf",
       list_of_significants, width = 10, height = 7)
ggsave("output/R/gene_content/enrichment_volcano.pdf",
       volcano, width = 15, height = 15)

# contingency_matrices["PAAR motif of membran proteins"]
# contingency_matrices["tail protein with lysin activity"]
