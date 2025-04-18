
classification %>%
  filter(Core == "yes") %>%
  count(Order, name = "contigs")

classified_taxa <- list()
core_taxonomy <- list()
# for (core_or_all in c("core", "all")) {
  for (tl in c("Order", "Family")) {
    tax_group <- paste0(tl, "_group")
    
    classif_filt <- classification
    # if (core_or_all == "core") {
    #   classif_filt <- classification %>%
    #     filter(Core == "yes") # %>%
    #     # mutate(Family_group = ifelse(Family_group == "ICTV-named_Caudoviricetes_family", Family, Family_group)) # This is breaking the logic a bit but it makes more sense here to list the actual family names
    # }
    
    core_taxonomy[[tl]] <- classification %>%
      filter(Core == "yes") %>%
      count(.data[[tl]], name = "contigs")
    
    classified_taxa[[tl]] <- classif_filt %>%
      select(starts_with(tl)) %>%
      filter(!str_detect(.data[[tax_group]], "Unclassified") & !str_detect(.data[[tax_group]], "unclassified")) %>%
      distinct() %>%
      group_by(.data[[tax_group]]) %>%
      summarise(taxa = n()) %>%
      rename("tax_group" = all_of(tax_group))
  }
# }

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
  arrange(desc(TPM)) %>%
  mutate(Class = factor(Class, levels = c("Caudoviricetes",
                                          "Faserviricetes|Huolimaviricetes|Malgrandaviricetes",
                                          "Vidaverviricetes|Leviviricetes",
                                          "Tokiviricetes",
                                          "Unclassified")))

pretty_pie$tibbles$load$Class <- phage_load$Class %>%
  pivot_longer(-Class) %>%
  group_by(Class) %>%
  mutate(load = mean(value)) %>%
  ungroup() %>%
  select(Class, load) %>%
  distinct() %>%
  arrange(desc(load)) %>%
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
  arrange(desc(TPM)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)

pretty_pie$tibbles$load$Order <- phage_load$Order_group %>%
  pivot_longer(-Order_group) %>%
  group_by(Order_group) %>%
  mutate(load = mean(value)) %>%
  ungroup() %>%
  select(Order_group, load) %>%
  distinct() %>%
  arrange(desc(load)) %>%
  mutate(Order_group = factor(Order_group, levels = c("Novel_Caudoviricetes_order", "Microviruses", "Novel_Tokiviricetes_order", "Unclassified"))) %>%
  rename(Order = Order_group)


pretty_pie$tibbles$n$Family <- classification %>%
  group_by(Family_group) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named_Caudoviricetes_family", "Novel_Tokiviricetes_family", "Unclassified_Microvirus", "Other_unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie$tibbles$TPM$Family <- phage_tpm$Family_group %>%
  pivot_longer(-Family_group) %>%
  group_by(Family_group) %>%
  mutate(TPM = mean(value)) %>%
  ungroup() %>%
  select(Family_group, TPM) %>%
  distinct() %>%
  arrange(desc(TPM)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named_Caudoviricetes_family", "Novel_Tokiviricetes_family", "Unclassified_Microvirus", "Other_unclassified"))) %>%
  rename(Family = Family_group)

pretty_pie$tibbles$load$Family <- phage_load$Family_group %>%
  pivot_longer(-Family_group) %>%
  group_by(Family_group) %>%
  mutate(load = mean(value)) %>%
  ungroup() %>%
  select(Family_group, load) %>%
  distinct() %>%
  arrange(desc(load)) %>%
  mutate(Family_group = factor(Family_group, levels = c("Novel_Caudoviricetes_family", "Microvirus_family", "ICTV-named_Caudoviricetes_family", "Novel_Tokiviricetes_family", "Unclassified_Microvirus", "Other_unclassified"))) %>%
  rename(Family = Family_group)

pretty_special_families <- list()

pretty_special_families$ictv$tibbles$n <- classification %>%
  filter(Family_group == "ICTV-named_Caudoviricetes_family") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))
pretty_special_families$micro$tibbles$n <- classification %>%
  filter(Family_group == "Microvirus_family") %>%
  group_by(Family) %>%
  count() %>%
  arrange(desc(n))

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
    labels <- NULL
    if (tl == "Class") {
      colors <- c("#8B4513", "#FFC300", "#D2691E", "#FFA07A", "#555555")
      labels <- levels(pretty_pie$tibbles[[thing]][[tl]][[tl]])
    }
    if (tl == "Order") {
      colors <- c("#FFDAB9", "#FFC300", "#FFA07A", "#555555")
      labels <- levels(pretty_pie$tibbles[[thing]][[tl]][[tl]]) %>%
        str_replace_all("_", " ")
    }
    if (tl == "Family") {
      colors <- c("#FFDAB9", "#FFC300", "#8B4513", "#FFA07A", "#666666", "#444444")
      labels <- levels(pretty_pie$tibbles[[thing]][[tl]][[tl]]) %>%
        str_replace_all("_", " ")
    }
    if (thing == "n") {
      counts_in_order <- pretty_pie$tibbles[[thing]][[tl]] %>%
        ungroup() %>%
        arrange(.data[[tl]]) %>%
        select(n) %>%
        unlist(use.names = FALSE)
      taxcount_in_order <- NULL
      if (!is.null(classified_taxa[[tl]])) {
        taxcount_in_order <- pretty_pie$tibbles[[thing]][[tl]] %>%
          ungroup() %>%
          arrange(.data[[tl]]) %>%
          left_join(., classified_taxa[[tl]], by = join_by(!!tl == tax_group)) %>%
          select(taxa) %>%
          unlist(use.names = FALSE)
      }
      labels <- paste0(labels, " (", counts_in_order, "/", taxcount_in_order, ")") %>%
        str_replace("/NA", "") %>%
        str_replace("/\\)", ")")
    }
  
    pretty_pie$plots[[thing]][[tl]] <- pretty_pie$tibbles[[thing]][[tl]] %>%
      ggplot(aes(x = "", y = .data[[thing]], fill = .data[[tl]])) +
      geom_bar(stat = "identity", color= "black") +
      # coord_polar(theta = "y", start = pi/2) +
      # coord_polar(theta = "y", start = 3/4 * pi) +
      # coord_polar(theta = "y", start = pi) +
      coord_polar(theta = "y", start = 7/4 * pi) +
      theme_void() +
      labs(fill = tl) +
      theme(legend.margin=margin(0,2,0,-20)) +
      scale_fill_manual(values = colors,
                        labels = labels) +
      ggtitle(thing)
    pretty_pie$plots[[thing]][[tl]] 
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
