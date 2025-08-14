library(tidyverse)

phold_per_cds_predictions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_per_cds_predictions_long_names.tsv.gz")

# The inputs are exported tables from Cytoscape.
microvirus_contigs <-  list()
microvirus_contigs[["0"]] <- read_lines("data/cytoscape_yFiles_Organic_layout_bin_0.topleft_blop")
for (subgraph in as.character(1:2)) { # <- mind the as.character()
  microvirus_contigs[[subgraph]] <- read.csv(paste0("output/vcontact3/bphage_vcontact3_b38_with_inphared/graph.bin_", subgraph, ".cyjs default node.csv")) %>%
    select(X_nx_name) %>%
    unlist(use.names = FALSE)
}
micros <- c()
for (subgraph in as.character(0:2)) {
  micros <- c(micros, microvirus_contigs[[subgraph]][str_detect(microvirus_contigs[[subgraph]], "NODE")])
    
}

min_evalue_repli <- phold_per_cds_predictions %>% 
  filter(contig_id %in% micros) %>% 
  filter(!is.na(evalue)) %>%
  group_by(contig_id) %>% 
  # filter("replication initiation protein" %in% product) %>%
  # filter("major head protein" %in% product) %>% 
  filter(product == "replication initiation protein") %>%  
  filter(evalue == min(evalue)) %>%
  ungroup()
min_evalue_head <- phold_per_cds_predictions %>% 
  filter(contig_id %in% micros) %>% 
  filter(!is.na(evalue)) %>%
  group_by(contig_id) %>% 
  # filter("replication initiation protein" %in% product) %>%
  # filter("major head protein" %in% product) %>%
  filter(product == "major head protein") %>% 
  filter(evalue == min(evalue)) %>%
  ungroup()

system("mkdir -p output/R/microviruses")

write_csv(min_evalue_repli, "output/R/microviruses/VP4_repli.csv")
write_csv(min_evalue_head, "output/R/microviruses/VP1_head.csv")

bind_rows(min_evalue_repli, min_evalue_head) %>%
  group_by(contig_id) %>% 
  filter("replication initiation protein" %in% product) %>%
  filter("major head protein" %in% product) %>%
  ungroup() %>%
  arrange(contig_id) %>% 
  write_csv("output/R/microviruses/micros_repli_and_head.csv")

write_lines(micros, "output/R/microviruses/microvirus_contigs")
