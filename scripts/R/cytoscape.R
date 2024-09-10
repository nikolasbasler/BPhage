library(tidyverse)

# classification <- readRDS("output/R/R_variables/classification.RDS")
# classification_gnmd <- readRDS("output/R/R_variables/classification_gnmd.RDS")

# The inputs are exported tables from Cytoscape.
graph_table <- list()
for (subgraph in as.character(0:3)) { # <- mind the as.character()
  graph_table[[subgraph]] <- read.csv(paste0("output/vcontact3/bphage_vcontact3_b38_with_inphared/graph.bin_", subgraph, ".cyjs default node.csv"))
}

present_in_all_countries <- read_lines("data/core_contigs.txt")


new_graph_table <- list()
for (subgraph in names(graph_table)) {
  new_graph_table[[subgraph]] <- graph_table[[subgraph]] %>%
    mutate(dataset = "Database") %>%
    mutate(dataset = ifelse(str_detect(X_nx_name, "NODE"), "BPhage", dataset)) %>%
    mutate(dataset = ifelse(X_nx_name %in% present_in_all_countries, "BPhage_core", dataset)) %>%
    mutate(dataset = ifelse(str_detect(X_nx_name, "Bonilla"), "Bonilla-Rosso", dataset)) %>%
    mutate(dataset = ifelse(str_detect(X_nx_name, "Deboutte"), "Deboutte", dataset)) %>%
    mutate(dataset = ifelse(str_detect(X_nx_name, "Busby"), "Busby", dataset)) %>%
    mutate(dataset = ifelse(!family == "", "Database_f", dataset)) %>%
    mutate(dataset = ifelse(!genus == "", "Database_g", dataset))
}

for (subgraph in names(new_graph_table)) {
  write_csv(new_graph_table[[subgraph]], paste0("output/vcontact3/bphage_vcontact3_b38_with_inphared/new_graph.bin_", subgraph, ".csv"))
}
