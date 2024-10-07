pretty_pie_tibble <- list()
pretty_pie_tibble_TPM <- list()

pretty_pie_tibble$Class <- classification %>% 
  group_by(Class) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes", 
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Unclassified",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes")))
pretty_pie_tibble_TPM$Class <- phage_tpm$Class %>%
  pivot_longer(-Class) %>%
  group_by(Class) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Class, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes", 
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Unclassified",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes")))

pretty_pie_tibble$Order <- classification %>%
  group_by(Order_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel", "Microviruses", "Unclassified"))) %>%
  rename(Order = Order_group)

pretty_pie_tibble_TPM$Order <- phage_tpm$Order_group %>%
  pivot_longer(-Order_group) %>%
  group_by(Order_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Order_group, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>% 
  mutate(Order_group = factor(Order_group, levels = c("Novel", "Microviruses", "Unclassified"))) %>%
  rename(Order = Order_group)


pretty_pie_tibble$Family <- classification %>% 
  group_by(Family_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel", "Microvirus", "ICTV-Named", "Unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie_tibble_TPM$Family <- phage_tpm$Family_group %>%
  pivot_longer(-Family_group) %>%
  group_by(Family_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family_group, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>% 
  mutate(Family_group = factor(Family_group, levels = c("Novel", "Microvirus", "ICTV-Named", "Unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie <- list()
pretty_pie_TPM <- list()
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

  pretty_pie_TPM[[tl]] <- pretty_pie_tibble_TPM[[tl]] %>%
    ggplot(aes(x = "", y = TPM, fill = .data[[tl]])) +
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
print(paste0(novel_families, " novel families"))
novel_orders <- classification %>% 
  select(Order) %>%
  filter(str_detect(Order, "novel")) %>%
  distinct() %>%
  nrow()
print(paste0(novel_orders, " novel orders"))
