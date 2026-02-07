library(tidyverse)
library(patchwork)

comps = c("comp_100", "comp_083", "comp_067", "comp_050")
names(comps) <- c("100", "083", "067", "050")

thresh = 0.9

all_predictions <- list()
for (comp in names(comps)) {
  bacphlip_predictions <- read.delim(paste0("output/lifestyle/tool_test/all_seqs_completeness_", comp, ".fasta.bacphlip")) %>%
    tibble() %>%
    rename(sample_name = X) %>%
    mutate(BACPHLIP = case_when(
      Virulent > Temperate & Virulent >= thresh ~ "Virulent",
      Virulent >= Temperate & Virulent < thresh ~ "Virulent (uncertain)",
      Temperate > Virulent & Temperate >= thresh ~ "Temperate",
      Temperate > Virulent & Temperate < thresh ~ "Temperate (uncertain)",
      
    )) %>%
    select(sample_name, BACPHLIP)
  
  replidec_predictions <- read.delim(paste0("output/lifestyle/tool_test/replidec_run_", comp, "_BC_predict.summary")) %>%
    tibble() %>%
    rename(Replidec = final_label) %>%
    select(sample_name, Replidec)
  
  phatyp_predictions <- read.delim(paste0("data/PhaTYP/phatyp_prediction_", comp,".tsv")) %>%
    tibble() %>%
    rename(
      sample_name = Accession,
      PhaTYP = TYPE
    ) %>%
    mutate(PhaTYP = case_when(
      PhaTYP == "virulent" &  PhaTYPScore >= thresh ~ "Virulent",
      PhaTYP == "virulent" &  PhaTYPScore < thresh ~ "Virulent (uncertain)",
      PhaTYP == "temperate" & PhaTYPScore >= thresh ~ "Temperate",
      PhaTYP == "temperate" & PhaTYPScore < thresh ~ "Temperate (uncertain)",
      PhaTYP == "-" ~ "no prediction"
    )) %>%
    select(sample_name, PhaTYP)
  
  all_predictions[[comps[comp]]] <- full_join(
    bacphlip_predictions, replidec_predictions, by = "sample_name"
    ) %>%
    full_join(., phatyp_predictions, by = "sample_name") %>%
    mutate(dataset = str_extract(sample_name, "^[A-Za-z]*"),
           ground_truth = case_when(
             dataset == "Engel" ~ "Virulent",
             dataset == "DSMZ" ~ "Virulent",
             dataset == "Bueren" ~ "Temperate",
           ),
           BACPHLIP = factor(BACPHLIP, levels = c("Chronic", "Temperate", "Temperate (uncertain)", "Virulent (uncertain)", "Virulent", "no prediction")),
           Replidec = factor(Replidec, levels = c("Chronic", "Temperate", "Temperate (uncertain)", "Virulent (uncertain)", "Virulent", "no prediction")),
           PhaTYP = factor(PhaTYP, levels = c("Chronic", "Temperate", "Temperate (uncertain)", "Virulent (uncertain)", "Virulent", "no prediction")),
    )
}

color_vector <- c("Virulent" = "#1C3A3A",
                  "Virulent (uncertain)" = "#4CB3B3",
                  "Chronic" =  "#8B4513", 
                  "Temperate" ="#FFC300", 
                  "Temperate (uncertain)" = "#FFDAB9",
                  "no prediction" = "grey"
)

plots <- list()
plots2 <- list()
Engel_plots <- list()
for (compo in names(comps)) {
  for (gt in c("Virulent", "Temperate")) {
    plots[[paste0(gt, "_", compo)]] <- all_predictions[[comps[compo]]] %>%
      filter(ground_truth == gt) %>%
      select(sample_name, BACPHLIP, PhaTYP, Replidec) %>%
      pivot_longer(-sample_name, names_to = "method", values_to = "call") %>%
      group_by(method, call) %>%
      mutate(count = n()) %>%
      ungroup() %>%
      select(-sample_name) %>%
      distinct() %>%
      rename(Assignment = call) %>%
      complete(method, Assignment, fill = list(count = 0L)) %>%
      
      ggplot(aes(x = method, y = count, fill = Assignment)) +
      geom_col() +
      labs(y = paste0("completeness: ", as.integer(compo), "%"),
           x = gt
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(face = "bold")
      ) +
      scale_fill_manual(values = color_vector)
    
    plots2[[paste0(gt, "_", compo)]] <- plots[[paste0(gt, "_", compo)]] + 
      labs(x = paste0("completeness: ", as.integer(compo), "%"),
           y = gt
      )
    
  }
  
  Engel_plots[[compo]] <- all_predictions[[comps[compo]]] %>%
    filter(dataset == "Engel") %>%
    select(-c(dataset, ground_truth)) %>%
    pivot_longer(-sample_name, names_to = "method", values_to = "call") %>%
    group_by(method, call) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    select(-sample_name) %>%
    distinct() %>%
    ggplot(aes(x = method, y = count, fill = call)) +
    geom_col() +
    ggtitle(paste0("Engel's phages\ncompleteness: ", as.integer(compo), "%")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = color_vector) +
    scale_y_continuous(breaks=c(seq(2,10,2)))
  
}

comparison_wrap <- wrap_plots(plots, 
                              guides = "collect", axes = "collect", 
                              labels = "collect", ncol = 2
)

plots2 <- list(
  Virulent_100 = plots2$Virulent_100, Virulent_083 = plots2$Virulent_083, Virulent_067 = plots2$Virulent_067, Virulent_050 = plots2$Virulent_050,
  Temperate_100 = plots2$Temperate_100, Temperate_083 = plots2$Temperate_083, Temperate_067 = plots2$Temperate_067, Temperate_050 = plots2$Temperate_050
)

comparison_wrap2 <- wrap_plots(plots2, 
                               guides = "collect", axes = "collect", 
                               labels = "collect", ncol = 4
)

engel_wrap <- wrap_plots(Engel_plots,
                         guides = "collect", axes = "collect", 
                         labels = "collect", ncol = 1)

full_tibble <- all_predictions$comp_100 %>%
  relocate(BACPHLIP, Replidec, PhaTYP,
           .after = ground_truth) %>%
  rename(
    BACPHLIP_100 = BACPHLIP,
    Replidec_100 = Replidec,
    PhaTYP_100 = PhaTYP
  ) %>%
  left_join(., all_predictions$comp_083[c("sample_name", "BACPHLIP", "Replidec", "PhaTYP")], by = "sample_name") %>%
  rename(
    BACPHLIP_083 = BACPHLIP,
    Replidec_083 = Replidec,
    PhaTYP_083 = PhaTYP
    
  ) %>%
  left_join(., all_predictions$comp_067[c("sample_name", "BACPHLIP", "Replidec", "PhaTYP")], by = "sample_name") %>%
  rename(
    BACPHLIP_067 = BACPHLIP,
    Replidec_067 = Replidec,
    PhaTYP_067 = PhaTYP
  ) %>%
  left_join(., all_predictions$comp_050[c("sample_name", "BACPHLIP", "Replidec", "PhaTYP")], by = "sample_name") %>%
  rename(
    BACPHLIP_050 = BACPHLIP,
    Replidec_050 = Replidec,
    PhaTYP_050 = PhaTYP
  )

#####
# Save files

system("mkdir -p output/R/lifestyle/tool_test/")
ggsave("output/R/lifestyle/tool_test/tool_comparison.pdf", comparison_wrap,
       height = 8, width = 6)
ggsave("output/R/lifestyle/tool_test/tool_comparison_horizontal.pdf", comparison_wrap2,
       height = 5, width = 10)
ggsave("output/R/lifestyle/tool_test/engels_phages.pdf", engel_wrap,
       height = 8, width = 3.5)

write_delim(full_tibble, "output/R/lifestyle/tool_test/tool_comparison.tsv",
            delim = "\t")
