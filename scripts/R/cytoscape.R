library(tidyverse)
# To run this script: 
# vcontact3 generates several .cyjs files. Open them in Cytoscape according to 
# their instructions on https://bitbucket.org/MAVERICLab/vcontact3/src/master/ 
# and export the tables by clicking on "Export Table to file..." above the 
# tabular display. Save them into 
# output/vcontact3/bphage_vcontact3_b38_with_inphared/ with the suggested
# default name. These csv tables are the inputs for this script.
# This script will then add another column to the table wich can then be 
# imported into Cytoscape again by clickin on "Import Table from file...". Now
# the nodes in the network can be colored according to the dataset by clicking
# on the "Style" tab at the left hand side then on "Fill Color". As Column
# select dataset, Mapping Type "Discrete Mapping" and choose colors for each
# dataset. 

# The colors chosen here are:
# BPhage: #FEC44F
# Bonilla-Rosso: #DD3497
# Busby: #41AB5D
# Database: #BDBDBD
# Database_f: #737373 # Database genome with family classification
# Database_g: #252525 # Database genome with genus classification
# Deboutte: #0570B0

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
