library(tidyverse)
library(ggpubr)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")

# phold_all_cds_functions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_all_cds_functions_long_names.tsv.gz")
phold_per_cds_predictions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_per_cds_predictions_long_names.tsv.gz")

present_in_all_countries <- read_lines("data/core_contigs.txt")



# SOMEHOW MAKE IT SO THAT ONLY BOXES WITH SIGNIFICANT PWCs ARE DRAWN 
# I hard-coded this now. Automate it.

# test_metrics <- c("function.", "product", "transl_table")
# annotation_core_plots <- list()
# for (metric in test_metrics) {
#   annotation_core_plots[[metric]] <- phold_per_cds_predictions %>%
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

test_metrics <- c("function.", "product", "transl_table")
annotation_core_plots <- list()
for (metric in test_metrics) {
  plot_tbl <- phold_per_cds_predictions %>%
    mutate(transl_table = as.character(transl_table)) %>%
    filter(str_starts(contig_id, "NODE")) %>%
    group_by(contig_id, .data[[metric]]) %>%
    count() %>%
    ungroup() %>% 
    left_join(., classification, by = join_by(contig_id == contig)) %>%
    filter(!is.na(Core),
           !is.na(completeness)) %>% 
    filter(completeness == 100) %>% mutate(normalised_gene_count = n) %>%
    # group_by(Core) %>%
    # mutate(normalised_gene_count = n / (completeness/100)) %>%
    # mutate(normalised_gene_count = n / length_kb) %>% ## <- this used for now.
    # mutate(genes_per_genome_kb = n / predicted_genome_length) %>%
    # ungroup() %>%
    # select(all_of(metric), normalised_gene_count, genes_per_genome_kb, Core) %>%
    # select(all_of(metric), normalised_gene_count, Core) %>%
    # pivot_longer(-c(all_of(metric), Core)) %>%
    group_by(.data[[metric]]) %>%
    filter(any(Core == "yes") & any(Core == "no")) %>%
    ungroup() %>%
    filter(.data[[metric]] != "hypothetical protein")

   annotation_core_plots[[metric]] <- plot_tbl %>%
     ggplot(aes(x = .data[[metric]], y = normalised_gene_count, fill = Core)) +
     # ggplot(aes(x = all_of(metric), y = n)) +
     geom_boxplot(aes(x = .data[[metric]], y = normalised_gene_count)) +
    # facet_wrap(~name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = metric) +
     geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE, group.by = "x.var")
   
   if (metric == "product") {
     significant_products <- c("anti-repressor Ant", "baseplate hub", "baseplate spike",
                               "baseplate wedge subunit", "beta-propeller repeat protein",
                               "DefenseFinder protein", "DNA binding protein", "DNA polymerase",
                               "GTP-binding domain", "head decoration", 
                               "head maturation protease", "head morphogenesis",
                               "head scaffolding protein", "head-tail adaptor", "holin",
                               "homing endonuclease", "integrase", "lipoprotein", "major head protein",
                               "Mu Gam-like end protection", "nucleotide kinase",
                               "portal protein", "RecT-like ssDNA annealing protein",
                               "replication initiation protein", "RIIB lysis inhibitor",
                               "Rz-like spanin", "single strand DNA binding protein",
                               "structural protein", "tail assembly chaperone", 
                               "tail chaperone protein", "tail completion or Neck1 protein",
                               "tail fiber protein", "tail length tape measure protein",
                               "tail protein", "tail sheath", "tail terminator", 
                               "terminase large subunit", "terminase small subunit",
                               "transcriptional repressor", "virion structural protein")
     plot_tbl_filt <- plot_tbl %>%
       filter(product %in% significant_products)
     plot_tbl_filt %>% select(product) %>% distinct() %>% unlist(use.names = FALSE) %>% sort()

     annotation_core_plots$significant_products <- plot_tbl_filt %>%
       ggplot(aes(x = .data[[metric]], y = normalised_gene_count, fill = Core)) +
       # ggplot(aes(x = all_of(metric), y = n)) +
       geom_boxplot(aes(x = .data[[metric]], y = normalised_gene_count)) +
       # facet_wrap(~name) +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       labs(title = metric) +
       geom_pwc(method="wilcox.test", label="p.adj.signif", hide.ns = TRUE, group.by = "x.var")
       
   }
}



### Save files ####
system("mkdir -p output/R/gene_content")
for (plot in names(annotation_core_plots)) {
  wid <- 12
  hig <- 8
  if (plot == "product") {
    wid <- 50
    hig <- 10
  }
  ggsave(paste0("output/R/gene_content/gene_content.", plot, ".pdf"), 
         annotation_core_plots[[plot]], width = wid, height = hig, limitsize = FALSE)
}
