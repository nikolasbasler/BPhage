library(tidyverse)


phold_all_cds_functions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_all_cds_functions.tsv")
phold_per_cds_predictions <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_per_cds_predictions.tsv")

present_in_all_countries <- read_lines("data/core_contigs.txt")


# unfoldables <- c("MN855733.1_Deboutte_CDS_0017", "MN855734.1_Deboutte_CDS_0010",
#                  "MN855779.1_Deboutte_CDS_0027", "NODE_A1_PT_19413_sum_rec_d_CDS_0015",
#                  "NODE_A2_RO_17377_sum_rec_d_CDS_0023", "NODE_A3_BE_16556_sum_rec_d_CDS_0022",
#                  "NODE_A3_BE_16562_aut_mid_d_CDS_0015", "NODE_A3_BE_16562_sum_rec_d_CDS_0074",
#                  "NODE_A3_NL_19103_aut_rec_d_CDS_0069", "NODE_A5_BE_16557_aut_rec_d_CDS_0034",
#                  "NODE_A6_BE_16562_aut_ile_d_CDS_0027")
# 
# phold_per_cds_predictions %>%
#   filter(cds_id %in% unfoldables) %>% View()


