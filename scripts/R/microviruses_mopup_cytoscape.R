library(tidyverse)

micros_id30_cytoscape_table <- read.csv("output/bphage_micros_mopup/bphage_micros_id30ForCytoscape.csv default node.csv")
micros_id50_cytoscape_table <- read.csv("output/bphage_micros_mopup/bphage_micros_id50ForCytoscape.csv default node.csv")

added_micros_id30_cytoscape_table <- micros_id30_cytoscape_table %>%
  mutate(subgroup = str_split_i(name, "_", 1)) %>% 
  mutate(subgroup = ifelse(subgroup == "NODE", "Bphage", subgroup))

added_micros_id50_cytoscape_table <- micros_id50_cytoscape_table %>%
  mutate(subgroup = str_split_i(name, "_", 1)) %>% 
  mutate(subgroup = ifelse(subgroup == "NODE", "Bphage", subgroup))

write_csv(added_micros_id30_cytoscape_table, "output/bphage_micros_mopup/bphage_micros_id30_cytoscape_table.csv")
write_csv(added_micros_id50_cytoscape_table, "output/bphage_micros_mopup/bphage_micros_id50_cytoscape_table.csv")
