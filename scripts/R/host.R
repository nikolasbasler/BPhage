##  LOAD lines 1-77 from alphabeta script


iphop_prediction_genus <- read.csv("output/host_prediction/iphop_output_core/Host_prediction_to_genus_m75.csv")
conservative_core_genera <- c("g__Bifidobacterium","g__Bombilactobacillus","g__Lactobacillus","g__Snodgrassella","g__Gilliamella")

iphop_prediction_genus %>%
  filter(Confidence.score >= 90) %>% 
  group_by(Virus) %>%
  filter(Confidence.score == max(Confidence.score)) %>%
  filter(grepl(paste(conservative_core_genera, collapse="|"), Host.genus)) %>%
  View()

  iphop_prediction_genus %>%
    select(Virus) %>%
    distinct()
  
  classification %>%
    filter(Host_group == "core" ) %>%
    select(Order) %>%
    distinct() %>% nrow()
  
  