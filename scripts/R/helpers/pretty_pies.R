pretty_pie_tibble <- list()
pretty_pie_tibble_TPM <- list()

pretty_pie_tibble$Class <- classification %>% 
  group_by(Class) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes", 
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes",
                                          "Unclassified")))
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
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes",
                                          "Unclassified")))

pretty_pie_tibble$Order <- classification %>%
  group_by(Order_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)

pretty_pie_tibble_TPM$Order <- phage_tpm$Order_group %>%
  pivot_longer(-Order_group) %>%
  group_by(Order_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Order_group, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>% 
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)


pretty_pie_tibble$Family <- classification %>% 
  group_by(Family_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named", "Novel_Tokiviricetes_family", "Unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie_tibble_TPM$Family <- phage_tpm$Family_group %>%
  pivot_longer(-Family_group) %>%
  group_by(Family_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family_group, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>% 
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named", "Novel_Tokiviricetes_family", "Unclassified"))) %>%
  rename(Family = Family_group)

ICTV_families <- list()
Microvirus_families <- list()

ICTV_families$n <- classification %>%
  filter(Family_group == "ICTV-named") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))
Microvirus_families$n <- classification %>%
  filter(Family_group == "Microvirus_family") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))


ICTV_ab <- phage_ab$Family %>% 
  filter(!str_starts(Family, "Microvirus")) %>%
  filter(!str_starts(Family, "novel")) %>%
  filter(Family != "Unclassified")
ICTV_families$TPM <- calc_tpm(abtable = ICTV_ab,
                              level = "Family", 
                              lengths_df = phage_lengths$Family) %>%
  select_if(~ all(!is.na(.))) %>%
  pivot_longer(-Family) %>%
  group_by(Family) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family, TPM) %>%
  distinct() %>%
  arrange(desc(TPM))

microvirus_ab <- phage_ab$Family %>% 
  filter(str_starts(Family, "Microvirus"))
Microvirus_families$TPM <- calc_tpm(abtable = microvirus_ab, 
                                    level = "Family", 
                                    lengths_df = phage_lengths$Family) %>%
  select_if(~ all(!is.na(.)))  %>%
  pivot_longer(-Family) %>%
  group_by(Family) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) 


pretty_pie <- list()
pretty_pie_TPM <- list()
for (tl in names(pretty_pie_tibble)) {
  if (tl == "Class") {
    colors <- c("#8B4513", "#FFC300", "#D2691E", "#FFA07A", "#555555")
    labels <- c("Caudoviricetes", 
                "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                "Vidaverviricetes|Leviviricetes",
                "Tokiviricetes",
                "Unclassified")
  }
  if (tl == "Order") {
    colors <- c("#FFDAB9", "#FFC300", "#FFA07A", "#555555")
    labels <- c("Novel Caudoviricetes order", "Microviruses", "Novel Tokiviricetes order", "Unclassified")
  }
  if (tl == "Family") {
    colors <- c("#FFDAB9", "#FFC300", "#8B4513", "#FFA07A", "#555555")
    labels <- c("Novel Caudoviricetes family", "Microvirus family", "ICTV-named", "Novel Tokiviricetes family", "Unclassified")
  }
  pretty_pie[[tl]] <- pretty_pie_tibble[[tl]] %>%
    ggplot(aes(x = "", y = n, fill = .data[[tl]])) +
    geom_bar(stat = "identity", color= "black") +
    coord_polar("y") +
    theme_void() +
    labs(fill = tl) +
    theme(legend.margin=margin(0,2,0,-20)) +
    scale_fill_manual(values = colors,
                      labels = labels)

  pretty_pie_TPM[[tl]] <- pretty_pie_tibble_TPM[[tl]] %>%
    ggplot(aes(x = "", y = TPM, fill = .data[[tl]])) +
    geom_bar(stat = "identity", color= "black") +
    coord_polar("y") +
    theme_void() +
    labs(fill = tl) +
    theme(legend.margin=margin(0,2,0,-20)) +
    scale_fill_manual(values = colors)
}

# microvirus_colors <- c("#FFF5CC", "#FFEB99", "#FFE066", "#FFD633", "#FFCC00", "#E6B800",
#                        "#CCA300", "#B38F00", "#997A00", "#806600", "#665200")
# ictv_colors <- c("#E6C2AA", "#D9B094", "#CC9D7F", "#BF8B6A", "#B27854", "#A5653F", "#99532A", 
#                  "#8B4513", "#7E3E11", "#70360F", "#632F0D", "#55270B", "#471F09")

ictv_colors <- c("#A5653F", "#471F09", "#B27854", "#55270B", "#BF8B6A", "#632F0D", "#99532A", 
                 "#CC9D7F", "#70360F", "#D9B094", "#7E3E11", "#E6C2AA", "#8B4513")
microvirus_colors <- rev(c("#FFCC00", "#665200", "#FFD633", "#806600", "#FFE066", "#997A00",
                           "#E6B800", "#FFEB99", "#B38F00", "#FFF5CC", "#CCA300"))

family_bar <- list()
for (thing in names(ICTV_families)) {
  family_bar$ictv[[thing]] <- ICTV_families[[thing]] %>%
    mutate(Family = factor(Family, levels = rev(ICTV_families[[thing]]$Family))) %>%
    ggplot(aes(x = "", y = .data[[thing]], fill = Family)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = ictv_colors) +
    theme_void() +
    theme(legend.margin=margin(0,10,0,-30)) 
    
  
  family_bar$micro[[thing]] <- Microvirus_families[[thing]] %>%
    mutate(Family = factor(Family, levels = rev(Microvirus_families[[thing]]$Family))) %>%
    ggplot(aes(x = "", y = .data[[thing]], fill = Family)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = microvirus_colors) +
    theme_void() +
    theme(legend.margin=margin(0,10,0,-30)) 
}


#

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
