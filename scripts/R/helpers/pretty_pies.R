pretty_pie <- list()

pretty_pie$tibbles$n$Class <- classification %>% 
  group_by(Class) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes", 
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes",
                                          "Unclassified")))
pretty_pie$tibbles$TPM$Class <- phage_tpm$Class %>%
  pivot_longer(-Class) %>%
  group_by(Class) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Class, TPM) %>%
  distinct() %>%
  # arrange(TPM) %>%
  # mutate(Class = factor(Class, levels = rev(c("Caudoviricetes", 
  #                                         "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
  #                                         "Vidaverviricetes|Leviviricetes",
  #                                         "Tokiviricetes",
  #                                         "Unclassified"))))
  arrange(desc(TPM)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes",
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes",
                                          "Unclassified")))

pretty_pie$tibbles$n$Order <- classification %>%
  group_by(Order_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)

pretty_pie$tibbles$TPM$Order <- phage_tpm$Order_group %>%
  pivot_longer(-Order_group) %>%
  group_by(Order_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Order_group, TPM) %>%
  distinct() %>%
  # arrange(TPM) %>% 
  # mutate(Order_group = factor(Order_group, levels = rev(c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified")))) %>%
  arrange(desc(TPM)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)


pretty_pie$tibbles$n$Family <- classification %>% 
  group_by(Family_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named", "Novel_Tokiviricetes_family", "Unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie$tibbles$TPM$Family <- phage_tpm$Family_group %>%
  pivot_longer(-Family_group) %>%
  group_by(Family_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family_group, TPM) %>%
  distinct() %>%
  # arrange(TPM) %>% 
  # mutate(Family_group = factor(Family_group, levels = rev(c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named", "Novel_Tokiviricetes_family", "Unclassified")))) %>%
  arrange(desc(TPM)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named", "Novel_Tokiviricetes_family", "Unclassified"))) %>%
  rename(Family = Family_group)

pretty_special_families <- list()

pretty_special_families$ictv$tibbles$n <- classification %>%
  filter(Family_group == "ICTV-named") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))
pretty_special_families$micro$tibbles$n <- classification %>%
  filter(Family_group == "Microvirus_family") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))

# This would re-calculate the TPM only for the selected families:
#
# ICTV_ab <- phage_ab$Family %>%
#   filter(!str_starts(Family, "Microvirus")) %>%
#   filter(!str_starts(Family, "novel")) %>%
#   filter(Family != "Unclassified")
# pretty_special_families$ictv$tibbles$TPM <- calc_tpm(abtable = ICTV_ab,
#                               level = "Family",
#                               lengths_df = phage_lengths$Family) %>%
#   select_if(~ all(!is.na(.))) %>%
#   pivot_longer(-Family) %>%
#   group_by(Family) %>%
#   mutate(TPM = mean(value)) %>%
#   ungroup() %>%
#   select(Family, TPM) %>%
#   distinct() %>%
#   arrange(desc(TPM))
# 
# microvirus_ab <- phage_ab$Family %>% 
#   filter(str_starts(Family, "Microvirus"))
# pretty_special_families$micro$tibbles$TPM <- calc_tpm(abtable = microvirus_ab, 
#                                     level = "Family", 
#                                     lengths_df = phage_lengths$Family) %>%
#   select_if(~ all(!is.na(.)))  %>%
#   pivot_longer(-Family) %>%
#   group_by(Family) %>%
#   mutate(TPM = mean(value)) %>%
#   ungroup() %>%
#   select(Family, TPM) %>%
#   distinct() %>%
#   arrange(desc(TPM)) 

pretty_special_families$ictv$tibbles$TPM <- phage_tpm$Family %>%
  filter(!str_starts(Family, "Microvirus")) %>%
  filter(!str_starts(Family, "novel")) %>%
  filter(Family != "Unclassified") %>%
  pivot_longer(-Family) %>%
  group_by(Family) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family, TPM) %>%
  distinct() %>%
  arrange(desc(TPM))

pretty_special_families$micro$tibbles$TPM <- phage_tpm$Family %>%
  filter(str_starts(Family, "Microvirus")) %>%
  pivot_longer(-Family) %>%
  group_by(Family) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family, TPM) %>%
  distinct() %>%
  arrange(desc(TPM))

for (thing in names(pretty_pie$tibbles)) {
  for (tl in names(pretty_pie$tibbles[[thing]])) {
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
    if (thing == "n") {
      counts_in_order <- pretty_pie$tibbles[[thing]][[tl]] %>%
        ungroup() %>%
        arrange(.data[[tl]]) %>%
        select(n) %>%
        unlist(use.names = FALSE)
      labels <- paste0(labels, " (", counts_in_order, ")")
    }
  
    pretty_pie$plots[[thing]][[tl]] <- pretty_pie$tibbles[[thing]][[tl]] %>%
      ggplot(aes(x = "", y = .data[[thing]], fill = .data[[tl]])) +
      geom_bar(stat = "identity", color= "black") +
      coord_polar(theta = "y", start = pi/2) +
      theme_void() +
      labs(fill = tl) +
      theme(legend.margin=margin(0,2,0,-20)) +
      scale_fill_manual(values = colors,
                        labels = labels) +
      ggtitle(thing)

  }
}

for (family_group in names(pretty_special_families)) {
  if (family_group == "ictv") {
    colors <- c("#A5653F", "#471F09", "#B27854", "#55270B", "#BF8B6A", "#632F0D", "#99532A", 
                "#CC9D7F", "#70360F", "#D9B094", "#7E3E11", "#E6C2AA", "#8B4513")
  }
  if (family_group == "micro") {
    colors <- rev(c("#FFCC00", "#665200", "#FFD633", "#806600", "#FFE066", "#997A00",
                    "#E6B800", "#FFEB99", "#B38F00", "#FFF5CC", "#CCA300"))
  }
  
  for (thing in names(pretty_special_families[[family_group]]$tibbles)) {
    pretty_special_families[[family_group]]$plots[[thing]] <- pretty_special_families[[family_group]]$tibbles[[thing]] %>%
      mutate(Family = factor(Family, levels = rev(pretty_special_families[[family_group]]$tibbles[[thing]]$Family))) %>%
      ggplot(aes(x = "", y = .data[[thing]], fill = Family)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = colors) +
      theme_void() +
      theme(legend.margin=margin(0,10,0,-10),
            axis.title.x = element_text(margin = margin(t = -20, b = 5))) +
      labs(x = thing)
  }
}

# Patchwork plots
pretty_patches <- list()
for (tl in names(pretty_pie$plots$n)) {
  pretty_patches[[tl]] <- pretty_pie$plots$n[[tl]] + 
    theme(
      legend.position = "top",        
      legend.direction = "horizontal", 
      legend.margin = margin(15, 0, 0, 250), 
      legend.box.margin = margin(-10, 0, 0, 0) 
    ) + 
    pretty_pie$plots$TPM[[tl]] + 
    theme(legend.position="none") 
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
