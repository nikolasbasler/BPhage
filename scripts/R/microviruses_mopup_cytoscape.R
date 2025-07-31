library(tidyverse)
# To run this script: Open the MOP-UP output file
# output/bphage_micros_mopup/bphage_micros_id30ForCytoscape.csv in Cytoscape,
# according their instructions on https://github.com/martinez-zacharya/MOP-UP
# and export the table by clicking on "Export Table to file..." above the 
# tabular display. Save it into output/bphage_micros_mopup/ with the suggested
# default name. That csv table is the inputs for this script.
# This script will then add another column to the table which can then be 
# imported into Cytoscape again by clickin on "Import Table from file...". Now
# the nodes in the network can be colored according to the microvirus family by 
# clicking on the "Style" tab at the left hand side then on "Fill Color". As 
# Column select subgroup, Mapping Type "Discrete Mapping" and choose colors for 
# each family.

# The colors chosen here are (families defined by Kirchberger et al 2022 are the
# same as in their paper):
# BPhage: #000000
# Family1: #646464
# Family10: #6F3984
# Family11: #9E37FF
# Family12: #2DB21B
# Family13: #FF67CC
# Family14: #B7A879
# Family15: #FEFC81
# Family16: #FFFFCC
# Family17: #855E9C
# Family18: #FFEF00
# Family19: #573D91
# Family2: #CCCD86
# Family3: #093190
# Family4: #FFD700
# Family4b: #FFD700
# Family5: #53BC7A
# Family5b: #53BC7A
# Family6: #FFFFCC
# Family7: #7C5835
# Family8: #855F9C
# Family9: #5D8AA8
# Obscuriviridae: #C8C8C8
# Singleton: #88CFF9


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
