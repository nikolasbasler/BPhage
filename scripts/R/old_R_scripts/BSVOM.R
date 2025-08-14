library(tidyverse)
classification <- readRDS("output/R/R_variables/classification.RDS")

pretty_pie_tibble <- list()
pretty_pie_tibble$Class <- classification %>% 
  group_by(Class) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes", 
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Unclassified",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes")))

pretty_pie_tibble$Order <- classification %>% 
  mutate(Order = ifelse(str_detect(Order, "novel_order"), "Novel", Order)) %>%
  group_by(Order) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Order = factor(Order, levels = c("Novel", "Microviruses", "Unclassified")))


pretty_pie_tibble$Family <- classification %>% 
  mutate(Family = ifelse(str_detect(Family, "novel_family"), "Novel", Family)) %>% 
  mutate(Family = ifelse(!Family == "Novel" & !str_detect(Family, "Micro") & Family !="Unclassified", "ICTV-Named", Family)) %>%
  mutate(Family = ifelse(str_detect(Family, "Micro"), "Microvirus", Family)) %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Family = factor(Family, levels = c("Novel", "Microvirus", "ICTV-Named", "Unclassified")))


pretty_pie <- list()
for (tl in names(pretty_pie_tibble)) {
  if (tl == "Class") {
   colors <- c("#8B4513", "#FFC300", "#555555", "#D2691E", "#FFA07A")
   }
  if (tl == "Order") {
    colors <- c("#FFDAB9", "#FFC300", "#555555")
  }
  if (tl == "Family") {
    colors <- c("#FFDAB9", "#FFC300", "#8B4513", "#555555")
  }
  pretty_pie[[tl]] <- pretty_pie_tibble[[tl]] %>%
    ggplot(aes(x = "", y = n, fill = .data[[tl]])) +
    geom_bar(stat = "identity", color= "black") +
    coord_polar("y") +
    theme_void() +
    labs(fill = tl) +
    theme(legend.margin=margin(0,2,0,-20)) +
    scale_fill_manual(values = colors)
}

novel_families <- classification %>% 
  select(Family) %>%
  filter(str_detect(Family, "novel")) %>%
  distinct() %>%
  nrow()
novel_orders <- classification %>% 
  select(Order) %>%
  filter(str_detect(Order, "novel")) %>%
  distinct() %>%
  nrow()


system("mkdir -p output/R/BSVOM_pretty_pies")
write_lines(novel_families, "output/R/BSVOM_pretty_pies/novel_families.txt")
write_lines(novel_orders, "output/R/BSVOM_pretty_pies/novel_orders.txt")
for (tl in names(pretty_pie_tibble)) {
  wid <- 4.5
  hei <- 4.5
  if (tl == "Class") {
    wid <- 6.6
    hei <- 4.5
  }
  write_csv(pretty_pie_tibble[[tl]],
            paste0("output/R/BSVOM_pretty_pies/pretty_pie.", tl,".csv"))
  ggsave(paste0("output/R/BSVOM_pretty_pies/pretty_pie.", tl,".pdf"),
         pretty_pie[[tl]], width = wid, height = hei)
}
 

