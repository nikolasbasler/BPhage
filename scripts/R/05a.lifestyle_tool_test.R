library(tidyverse)
library(patchwork)

comps = c("comp_100", "comp_083", "comp_067", "comp_050")
names(comps) <- c("100", "083", "067", "050")

bacphlip_threshold = 0.9

both_predictions <- list()
for (comp in names(comps)) {
  bacphlip_predictions <- read.delim(paste0("output/lifestyle/tool_test/all_seqs_completeness_", comp, ".fasta.bacphlip")) %>%
    tibble() %>%
    rename(sample_name = X) %>%
    mutate(
      bacphlip_lenient = ifelse(Virulent > Temperate, "Virulent", "Temperate"),
      bacphlip_stringent = case_when(
        Virulent >= bacphlip_threshold ~ "Virulent",
        Temperate >= bacphlip_threshold ~ "Temperate",
        .default = "uncertain"
      )
    ) %>%
    select(sample_name, bacphlip_lenient, bacphlip_stringent)
  
  replidec_predictions <- read.delim(paste0("output/lifestyle/tool_test/replidec_run_", comp, "_BC_predict.summary")) %>%
    tibble() %>%
    rename(replidec_final_label = final_label) %>%
    select(sample_name, replidec_final_label)
  
  both_predictions[[comps[comp]]] <- full_join(
    bacphlip_predictions, 
    replidec_predictions, 
    by = "sample_name") %>%
    mutate(dataset = str_extract(sample_name, "^[A-Za-z]*"),
           ground_truth = case_when(
             dataset == "Engel" ~ "Virulent",
             dataset == "DMSZ" ~ "Virulent",
             dataset == "Bueren" ~ "Temperate"
           )
    )
}

color_vector <- c("Virulent" = "#1C3A3A", 
                  "Chronic" =  "#8B4513", 
                  "Temperate" ="#FFC300", 
                  "uncertain" = "#777777"
                  )


plots <- list()
Engel_plots <- list()
for (compo in names(comps)) {
  for (gt in c("Virulent", "Temperate")) {
    plots[[paste0(gt, "_", compo)]] <- both_predictions[[comps[compo]]] %>%
      filter(ground_truth == gt) %>%
      select(sample_name, bacphlip_lenient, bacphlip_stringent, replidec_final_label) %>%
      pivot_longer(-sample_name, names_to = "method", values_to = "call") %>%
      group_by(method, call) %>%
      mutate(count = n()) %>%
      ungroup() %>%
      select(-sample_name) %>%
      distinct() %>%
      mutate(method = case_when(
        method == "bacphlip_lenient" ~ "BACPHLIP\n(lentient)",
        method == "bacphlip_stringent" ~ "BACPHLIP\n(stringent)",
        method == "replidec_final_label" ~ "Replidec"
        )) %>%
      rename(Assignment = call) %>%
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
  }
  
  Engel_plots[[compo]] <- both_predictions[[comps[compo]]] %>%
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

engel_wrap <- wrap_plots(Engel_plots,
           guides = "collect", axes = "collect", 
           labels = "collect", ncol = 1)

full_tibble <- both_predictions$comp_100 %>%
  relocate(bacphlip_lenient, bacphlip_stringent, replidec_final_label,
           .after = ground_truth) %>%
  rename(
    bacphlip_lenient_100 = bacphlip_lenient,
    bacphlip_stringent_100 = bacphlip_stringent,
    replidec_final_label_100 = replidec_final_label
    ) %>%
  left_join(., both_predictions$comp_083[c("sample_name", "bacphlip_lenient", "bacphlip_stringent", "replidec_final_label")], by = "sample_name") %>%
  rename(
    bacphlip_lenient_83 = bacphlip_lenient,
    bacphlip_stringent_83 = bacphlip_stringent,
    replidec_final_label_83 = replidec_final_label
  ) %>%
  left_join(., both_predictions$comp_067[c("sample_name", "bacphlip_lenient", "bacphlip_stringent", "replidec_final_label")], by = "sample_name") %>%
  rename(
    bacphlip_lenient_67 = bacphlip_lenient,
    bacphlip_stringent_67 = bacphlip_stringent,
    replidec_final_label_67 = replidec_final_label
  ) %>%
  left_join(., both_predictions$comp_050[c("sample_name", "bacphlip_lenient", "bacphlip_stringent", "replidec_final_label")], by = "sample_name") %>%
  rename(
    bacphlip_lenient_50 = bacphlip_lenient,
    bacphlip_stringent_50 = bacphlip_stringent,
    replidec_final_label_50 = replidec_final_label
  )

#####
# Save files

system("mkdir -p output/R/lifestyle/tool_test/")
ggsave("output/R/lifestyle/tool_test/tool_comparison.pdf", comparison_wrap,
       height = 8, width = 6)
ggsave("output/R/lifestyle/tool_test/engels_phages.pdf", engel_wrap,
       height = 8, width = 3.5)

write_delim(full_tibble, "output/R/lifestyle/tool_test/tool_comparison.tsv",
            delim = "\t")
